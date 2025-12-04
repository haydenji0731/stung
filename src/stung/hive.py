import numpy as np
from stung import bumble
from collections import defaultdict
import os

def dissect(
    mat : np.ndarray,
    n : int,
    out_dir : str,
    stung_blck : bumble.Block,
):
    """
    annotates various genomic rearrangements events including 
    insertions, duplications, non-colinear matches, and long
    stretches of inverse diagonal matches
    args : 
        mat (np.ndarray) : 2d match matrix
        n (int) : dimension
        out_dir (str) : output directory
        stung_blck (bumble.Block) : a stung (= breakpoint in between colienar blocks) instance
    returns :
        None
    raises :
        None
    """
    y_start = stung_blck.start.y
    y_end = stung_blck.end.y
    x_start = stung_blck.start.x
    x_end = stung_blck.end.x

    # annotate same strand matches
    matches = set() # (x, y) coordinates

    y = y_start
    x = x_start
    while y < y_end and x < x_end:
        if mat[y][x] == 1:
            matches.add((x, y))
        else:
            break
        x += 1
        y += 1
    
    y = y_end
    x = x_end
    while y >= y_start and x >= x_start:
        if mat[y][x] == 1:
            matches.add((x, y))
        else:
            break
        x -= 1
        y -= 1

    matches = len(matches)
    matches.sort(key = lambda xy : xy[0])

    # annotate gaps (= where insertions occur)
    x_gaps = []
    y_gaps = []

    prev_x, prev_y = matches[0]
    for i in range(1, len(matches), 1):
        x, y = matches[i]
        if x - prev_x > 1:
            # gap_start, gap_length
            x_gaps.append((prev_x + 1, x - prev_x - 1))
        if y - prev_y > 1:
            y_gaps.append((prev_y + 1, y - prev_y - 1))
        prev_x = x
        prev_y = y
    
    # annotate inverse diagonal matches (continuous)
    inv_matches, _ = find_inverse_diagonals(
        mat = mat,
        x_start = x_start,
        x_end = x_end,
        y_start = y_start,
        y_end = y_end,
        n = n
    )

    # annotate duplications (row and col separately)
    # also annotates additional insertions based on row and col sums
    row_dups, row_matches, col_dups, col_matches, x_ins, y_ins = find_duplications(
        mat = mat,
        diag_matches = matches,
        x_start = x_start,
        x_end = x_end,
        y_start = y_start,
        y_end = y_end
    )

    # overlay gaps, opposite strand matches, and duplications
    row_anns = [[] for _ in range(n)]
    col_anns = [[] for _ in range(n)]

    for start_i, l in x_gaps:
        for i in range(l):
            col_anns[start_i + i].append('INS')
    
    for start_i, l in y_gaps:
        for i in range(l):
            row_anns[start_i + i].append('INS')
    
    for x, y in inv_matches:
        col_anns[x].append(f'INV_BLCK:{y}')
        row_anns[y].append(f'INV_BLCK:{x}')
    
    for y in row_dups:
        for x, is_inverse in row_dups[y]:
            if is_inverse:
                row_anns[y].append(f'INV_DUP:{x}')
            else:
                row_anns[y].append(f'DUP:{x}')
    
    for y in row_matches:
        for x, is_inverse in row_matches[y]:
            if is_inverse:
                row_anns[y].append(f'INV_MATCH:{x}')
            else:
                row_anns[y].append(f'MATCH:{x}')
    
    for x in col_dups:
        for y, is_inverse in col_dups[x]:
            if is_inverse:
                col_anns[x].append(f'INV_DUP:{y}')
            else:
                col_anns[x].append(f'DUP:{y}')
    
    for x in col_matches:
        for y, is_inverse in col_matches[x]:
            if is_inverse:
                col_anns[x].append(f'INV_MATCH:{y}')
            else:
                col_anns[x].append(f'MATCH:{y}')
    
    for y, _ in y_ins:
        if 'INS' not in row_anns[y]:
            row_anns[y].append('INS')
    for x, _ in x_ins:
        if 'INS' not in col_anns[x]:
            row_anns[x].append('INS')
        
    # write results to file(s)
    _, _ = write_ann2file(
        out_dir = out_dir,
        n = n,
        row_anns = row_anns,
        col_anns = col_anns
    )

def write_ann2file(
    out_dir : str,
    n : int,
    row_anns : list[list[str]],
    col_anns : list[list[str]]
) -> tuple[str, str]:
    """
    writes genomic event annotations as TSV files
    args :
        out_dir (str) : output directory
        n (int) : dimension
        row_anns (list(list(str))) : list of row event annotations; each row is associated with a list of events
        col_anns (list(list(str))) : list of col event annotations; each col is associated with a list of events
    returns :
        row_file_path (str) : path to the file where row event annotations are saved
        col_file_path (str) : path to the file where col event annotations are saved
    raises :
        None
    """
    row_file_path = os.path.join(out_dir, 'row_events.ann')
    with open(row_file_path, 'w') as fh:
        for i in range(n):
            fh.write(f'{i}\t{";".join(row_anns[i])}\n')
    
    col_file_path = os.path.join(out_dir, 'col_events.ann')
    with open(col_file_path, 'w') as fh:
        for i in range(n):
            fh.write(f'{i}\t{";".join(col_anns[i])}\n')
    return row_file_path, col_file_path

def find_inverse_diagonals(
    mat : np.ndarray,
    x_start : int,
    x_end : int,
    y_start : int,
    y_end : int,
    n : int,
    min_l : int = 2
):
    """
    finds diagonal stretches of consecutive inverse matches
    args :
        mat (np.ndarray) : 2d match matrix
        x_start (int) : start coordinate on x-axis
        x_end (int) : end coordinate on x-axis
        y_start (int) : start coordinate on y-axis
        y_end (int) : end coordinate on y-axis
        n (int) : dimension
        min_l : minimum diagonal length to be considered a diagonal match
    returns :
        inv_matches (list(tuple(int, int))) : list of points that comprises inverse diagonal matches
        inv_diag_blocks (list(bumble.Block)) : list of inverse diagonal blocks
    raises :
        None
    """
    diag_mat = np.zeros((n, n))

    for i in range(y_end - 1, y_start - 1, -1):
        for j in range(x_start, x_end, 1):
            y, x = i, j
            ctr = 0
            while y >= y_start and x < x_end:
                if mat[y][x] == 2: # only look for opposite strand matches
                    ctr += 1
                else:
                    break
                y -= 1
                x += 1
            diag_mat[i][j] = ctr if ctr >= min_l else 0
    
    inv_matches = []
    inv_diag_blocks = []
    for i in range(y_end - 1, y_start - 1, -1):
        for j in range(x_start, x_end, 1):
            if diag_mat[i][j] > 0:
                l = int(diag_mat[i][j])
                for k in range(l):
                    inv_matches.append((j + k, i - k))
                    diag_mat[i - k][j + k] = 0
                inv_diag_blocks.append(bumble.Block(
                    start = bumble.Point(x = j, y = i),
                    end = bumble.Point(x = j + l, y = i - l)
                ))

    return inv_matches, inv_diag_blocks
        
def find_duplications(
    mat : np.ndarray,
    diag_matches : list[tuple[int, int]],
    x_start : int,
    x_end : int,
    y_start : int,
    y_end : int
):
    """
    classifies same strand and opposite strand matches as either
    diagonal (~colinear) or non-diagonal (~non-colinear) or duplications
    if a match is NOT part of a larger diagonal stretch, then it is either a non-diagonal
    match or a duplication depending whether a match already exists on some diagonal
    args :
        mat (np.ndarray) : 2d matrix
        diag_matches (list(tuple(int, int))) : list of diagonal matches
        x_start : start coordinate on the x-axis
        x_end : end coordinate on the x-axis
        y_start : start coordinate on the y-axis
        y_end : end coordinate on the y-axis
    returns :
        row_dups
        row_matches
        col_dups
        col_matches
        x_ins
        y_ins
    raises :
        None
    """
    row_dups = defaultdict(list) # h1 dups found on h0
    row_matches = defaultdict(list) # h1 off-diag matches found on h0
    col_dups = defaultdict(list) # h0 dups found on h1
    col_matches = defaultdict(list) # h0 off-diag matches found on h1

    y_ins = []
    x_ins = []
    
    for i in range(y_start, y_end, 1):
        curr_matches = np.where(mat[i] > 0)[0]
        if len(curr_matches) == 0:
            y_ins.append((i, 1))
            continue
        
        diag_match_found = False
        for j in curr_matches:
            if (j, i) in diag_matches:
                diag_match_found = True
                break

        for j in curr_matches:
            if (j, i) in diag_matches:
                continue
            if diag_match_found:
                row_dups[i].append((int(j), bool(mat[i][j] == 2)))
            else:
                row_matches[i].append((int(j), bool(mat[i][j] == 2)))
    
    for j in range(x_start, x_end, 1):
        curr_matches = np.where(mat[:, j] > 0)[0]
        if len(curr_matches) == 0:
            x_ins.append((j, 1))
            continue

        diag_match_found = False
        for i in curr_matches:
            if (j, i) in diag_matches:
                diag_match_found = True
                break
        
        for i in curr_matches:
            if (j, i) in diag_matches:
                continue
            if diag_match_found:
                col_dups[j].append((int(i), bool(mat[i][j] == 2)))
            else:
                col_matches[j].append((int(i), bool(mat[i][j] == 2)))

    return row_dups, row_matches, col_dups, col_matches, x_ins, y_ins