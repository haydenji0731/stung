import numpy as np
from stung import bumble

def annotate(
    mat,
    stung_blck : bumble.Block
):
    """
    TODO
    """
    y_start = stung_blck.start.y
    y_end = stung_blck.end.y
    x_start = stung_blck.start.x
    x_end = stung_blck.end.x

    # first annotate fwd matches
    matches = [] # (x, y) coordinates

    y = y_start
    x = x_start
    while y < y_end and x < x_end:
        if mat[y][x] == 1:
            matches.append((x, y))
        else:
            break
        x += 1
        y += 1
    
    y = y_end
    x = x_end
    while y >= y_start and x >= x_start:
        if mat[y][x] == 1:
            matches.append((x, y))
        else:
            break
        x -= 1
        y -= 1

    matches.sort(key = lambda xy : xy[0])

    # annotate gaps

    x_gaps = []
    y_gaps = []

    prev_x, prev_y = matches[0]
    for i in range(1, len(matches), 1):
        x, y = matches[i]
        if x - prev_x > 1:
            # gap_start, gap_length
            x_gaps.append((prev_x, x - prev_x - 1))
        if y - prev_y > 1:
            y_gaps.append((prev_y, y - prev_y - 1))
    
def find_inverse_diagonals(
    mat,
    x_start : int,
    x_end : int,
    y_start : int,
    y_end : int,
    n : int,
    min_l : int = 2
):
    """
    TODO
    """
    diag_mat = np.zeros((n, n))

    for i in range(y_end - 1, y_start - 1, -1):
        for j in range(x_start, x_end, 1):
            y, x = i, j
            ctr = 0
            while y >= y_start and x < x_end:
                if mat[y][x] == 2: # only look for rev matches
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
                inv_diag_blocks.append(
                    
                )

    return
        

