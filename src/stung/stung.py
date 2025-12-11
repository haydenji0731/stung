import numpy as np
from stung import bumble
from stung import hive
import os
import copy

def compute_mlen(
    mat : np.ndarray,
    n : int
) -> np.ndarray:
    """
    computes the maximum 1-filled positive diagonal length from each cell
    args :
        mat (np.ndarray) : 2d match matrix
        n (int) : match matrix dimension
    returns :
        diag_mlen_mat (np.ndarray) : a 2d matrix storing diagonal match lengths
    raises :
        None
    """
    diag_mlen_mat = np.zeros((n, n), dtype=int)
    for i in range(n): # row index
        for j in range(n): # col index
            mlen = 0
            y, x = i, j
            while x < n and y < n:
                if mat[y][x] != 1: # TODO : is this correct?
                    break
                mlen += 1
                x += 1
                y += 1
            diag_mlen_mat[i][j] = mlen
    return diag_mlen_mat

def find_colinear_blocks(
    mlen_mat : np.ndarray,
    n : int,
    min_dlen : int,
    verbose : bool = False
) -> list[bumble.Block]:
    """
    finds colinear blocks using diagonal match length matrix
    only blocks with diagonal lenght > min_dlen are returned
    args :
        mlen_mat (np.ndarray) : deep copy of the original mlen_mat; changes not propagate
        n (int) :  mlen_mat length
        min_dlen (int) : minimum diagonal length to qualify as a colinear block
        verbose (bool) : print informative messages
    returns :
        blcks (list(bumbl.Block)) : list of colinear blocks
    raises :
        None
    """

    blcks = []
    for i in range(n): # row index
        for j in range(n): # col index
            start = None
            end = None
            dlen = 0
            y, x = i, j
            while x < n and y < n:
                if mlen_mat[y][x] > 0:
                    if not start:
                        start = bumble.Point(x = x, y = y)
                    mlen_mat[y][x] = 0
                    dlen += 1
                else:
                    if start and dlen > 1:
                        end = bumble.Point(x = x, y = y)
                    break
                x += 1
                y += 1
            if start and not end and dlen > 1:
                end = bumble.Point(x = x, y = y)
            
            if not start or not end:
                continue
        
            if dlen < min_dlen:
                continue

            blcks.append(
                bumble.Block(
                    start = start, 
                    end = end)
            )
    if verbose:
        print(f'{len(blcks)} colinear blocks of length > {min_dlen} detected')
    
    return blcks

def collapse_blocks(
    blocks = list[bumble.Block]
) -> tuple[bool, list[bumble.Block]]:
    """
    collapses neighboring blocks if two blocks' area of intersection cover
    at least 50% of one of them
    args :
        blocks (list(bumble.Block)) : input list of colinear blocks to collapse
    returns :
        ctr (int) : number of blocks post-collapse
        collapsed (list(bumble)) : processed list of colinear blocks post-collapse
    raises :
        None
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
            anchor = bumble.Block(
                anchor.start, 
                bumble.Point(
                    x = max(anchor.end.x, b.end.x),
                    y = max(anchor.end.y, b.end.y)
                )
            )
        else:
            anchor = b
            collapsed.append(b)
        
    return ctr, collapsed
            
def get_stungs(
    colin_blocks : list[bumble.Block],
    verbose : bool = False
) -> list[tuple[bumble.Point, bumble.Point]]:
    """
    get stungs given a sorted list of colinear blocks
    a stung is the region in between two consecutive colinear blocks
    a stung is represented a tuple of bumble.Point objects
    args :
        colin_blocks (list(bumble.Block)) : input sorted list of colinear blocks
        verbose (bool; default=False) : print informative messages
    returns :
        stungs (list(tuple(bumble.Point, bumble.Point))) : list of stungs whose coordinates are lcro
    raises :
        None
    """
    stungs = []
    for i in range(0, len(colin_blocks) - 1, 1):

        curr_blck = colin_blocks[i]
        next_blck = colin_blocks[i + 1]

        xst = curr_blck.end.x # inclusive
        xen = next_blck.start.x # exclusive
        yst = curr_blck.end.y # inclusive
        yen = next_blck.start.y # exclusive
        
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
    bmbl : bumble.Bumble,
    mat : np.ndarray,
    n : int,
    out_dir : str,
    verbose : bool = False
):
    """
    give a Bumble instance and a match matrix,
    compute colinear blocks and stungs before and after colinear block extensions
    based on pairwise sequence alignments
    genomic events in each stung are concurrently annotated and saved (see hive module for more info)
    args :
        bmbl (bumble.Bumble) : a Bumble instance storing initial settings for a stung analysis
        mat (np.ndarray) : 2d match matrix
        n (int) : match matrix length
        out_dir (str) : output directory
        verbose (bool; default=False) : print informative messages
    returns :
        init_mat (np.ndarray) : safe copy of the input match matrix
        stungs (list(tuple(bumble.Point, bumble.Point))) : list of stungs as tuples
        stung_blcks (list(bumble.Block)) : list of stungs as bumble.Block instances
    raises :
        ValueError
    """
    if len(mat) != n:
        raise ValueError(f'invalid dimension {n}')
    
    init_mat = copy.deepcopy(mat) # keep a copy of the initial match matrix

    diag_mlen_mat = compute_mlen(mat = mat, n = n)
    
    colin_blocks = find_colinear_blocks(
        mlen_mat = diag_mlen_mat.copy(),
        n = n,
        min_dlen = bmbl.min_dlen,
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

    # perform pairwise nucleotide sequence alignments
    if verbose:
        print("executing pairwise alignments around stungs...")

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

    # update colin_blocks
    diag_mlen_mat_2 = compute_mlen(mat = mat, n = n)
    colin_blocks_2 = find_colinear_blocks(
        mlen_mat = diag_mlen_mat_2.copy(),
        min_dlen = bmbl.min_dlen,
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