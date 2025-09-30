import numpy as np
from stung import bumble
from stung import hive
import os
import copy

def compute_mlen(
    mat,
    n : int
):
    diag_mlen_mat = np.zeros((n, n), dtype=int)
    for i in range(n): # row index
        for j in range(n): # col index
            mlen = 0
            y, x = i, j
            while x < n and y < n:
                if mat[y][x] < 1:
                    break
                mlen += 1
                x += 1
                y += 1
            diag_mlen_mat[i][j] = mlen
    return diag_mlen_mat

def find_colinear_blocks(
    mlen_mat,
    n : int,
    min_blen : int = 10,
    verbose : bool = False
) -> list[bumble.Block]:
    """
    TODO
    """

    blcks = []
    for i in range(n): # row index
        for j in range(n): # col indx
            start = None
            end = None
            blen = 0
            y, x = i, j
            while x < n and y < n:
                if mlen_mat[y][x] > 0: # TODO: is this okay?
                    if not start:
                        start = bumble.Point(x = x, y = y)
                    mlen_mat[y][x] = 0
                    blen += 1
                else:
                    if start and blen > 1:
                        end = bumble.Point(x = x, y = y)
                    break
                x += 1
                y += 1
            if start and not end and blen > 1:
                end = bumble.Point(x = x, y = y)
            
            if not start or not end:
                continue
        
            if blen < min_blen:
                continue

            blcks.append(
                bumble.Block(
                    start = start, 
                    end = end)
            )
    if verbose:
        print(f'{len(blcks)} colinear blocks of length > {min_blen} detected')
    return blcks

def collapse_blocks(
    blocks = list[bumble.Block]
) -> tuple[bool, list[bumble.Block]]:
    """
    TODO
    """
    l = len(blocks)
    collapsed = []
    used = [0] * l
    anchor = blocks[0]
    used[0] = True
    collapsed.append(anchor)
    ctr = 0

    for i in range(l):
        b = blocks[i]
        if used[i]:
            continue
        if bumble.do_connect_greedy(anchor, b):
            ctr += 1
            anchor = bumble.Block(anchor.start, b.end)
        else:
            anchor = b
            collapsed.append(b)
        
    return ctr, collapsed
            
def get_stungs(
    colin_blocks : list[bumble.Block],
    verbose : bool = False
) -> list[tuple[bumble.Point, bumble.Point]]:
    """
    TODO
    """
    stungs = []
    for i in range(0, len(colin_blocks) - 1, 1):

        curr_blck = colin_blocks[i]
        next_blck = colin_blocks[i + 1]

        xst = curr_blck.end.x
        xen = next_blck.start.x
        yst = curr_blck.end.y
        yen = next_blck.start.y
        
        stungs.append((
            bumble.Point(
                x = xst,
                y = yst
            ),
            bumble.Point(
                x = xen,
                y = yen
            )
        ))
    if verbose:
        print(f'{len(stungs)} (st)ructurally (un)stable (g)enomic regions detected')
    
    return stungs

def buzz(
    bmbl,
    mat,
    n : int,
    out_dir : str,
    verbose : bool = False
):
    """
    TODO
    """
    if len(mat) != n:
        raise ValueError()
    
    init_mat = copy.deepcopy(mat)

    diag_mlen_mat = compute_mlen(mat = mat, n = n)
    
    colin_blocks = find_colinear_blocks(
        mlen_mat = diag_mlen_mat.copy(), # pass a shallow copy to prevent mods from propagating
        n = n,
        verbose = verbose
    )

    n_collapsed, colin_blocks = collapse_blocks(
        blocks = colin_blocks
    )
    if verbose:
        print(f'{n_collapsed} blocks were collapsed')

    bumble.write_blocks2file(colin_blocks, os.path.join(out_dir, 'pre_colin_blcks.tsv'))
    
    stungs = get_stungs(
        colin_blocks = colin_blocks,
        verbose = verbose
    )

    stung_blcks = []
    for i, stung in enumerate(stungs):
        xlim, ylim, wkdir = bmbl.extend_colinear_blocks(
            mat = mat,
            n = n,
            stung = stung,
            bp_i = i
        )
        stung_blck = bumble.Block(
            start = bumble.Point(
                x = xlim[0],
                y = ylim[0]
            ),
            end = bumble.Point(
                x = xlim[1],
                y = ylim[1]
            )
        )
        stung_blcks.append(stung_blck)

        hive.dissect(
            mat = mat,
            n = n,
            out_dir = wkdir,
            stung_blck = stung_blck,
        )
        # break # remove after testing

    # update colin_blocks
    diag_mlen_mat_2 = compute_mlen(mat = mat, n = n)
    colin_blocks_2 = find_colinear_blocks(
        mlen_mat = diag_mlen_mat_2.copy(),
        n = n,
        verbose = verbose
    )
    n_collapsed_2, colin_blocks_2 = collapse_blocks(
        blocks = colin_blocks_2
    )
    if verbose:
        print(f'{n_collapsed_2} blocks were collapsed (post-alignment)')

    bumble.write_blocks2file(colin_blocks_2, os.path.join(out_dir, 'post_colin_blcks.tsv'))
    
    return init_mat, stungs, stung_blcks
