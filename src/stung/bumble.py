import pyfastx
import stung.utils as utl
import os
import numpy as np

class Point:
    def __init__(
        self,
        x : int,
        y : int
    ):
        self.x = x
        self.y = y

    def __str__(
        self
    ):
        return f"{self.x},{self.y}"

class Block:
    def __init__(
        self,
        start : Point,
        end : Point
    ):
        self.start = start
        self.end = end
        self.dx = end.x - start.x
        self.dy = end.y - start.y
        self.diag_len = ((self.dx) ** 2 + (self.dy) ** 2) ** 0.5
        self.area = self.dx * self.dy
    
    def __str__(
        self
    ):
        return f'{self.start}\t{self.end}'

def do_connect_greedy(
    anchor : Block,
    other : Block,
    max_ratio : float = 0.5
):
    """
    """
    dx = anchor.end.x - other.start.x
    dy = anchor.end.y - other.start.y
    dx = max(0, dx)
    dy = max(0, dy)
    A = dx * dy
    return A / other.area > max_ratio

def write_blocks2file(
    blocks : list[Block],
    out_fp : str
):
    """
    """
    with open(out_fp, 'w') as fh:
        for b in blocks:
            fh.write(f'{str(b)}\n')

# TODO: finish implementation
class Gene:
    def __init__(
        self
    ):
        raise NotImplementedError

class Bumble:
    def __init__(
        self,
        genome_fp : str,
        ann_fps : list[str],
        temp_dir : str,
        out_dir : str,
        min_pident : float = 90.0,
        min_diag_len : int = 10,
        pad_len : int = 10
    ):
        """
        TODO
        """
        try:
            self.genome = pyfastx.Fasta(genome_fp)
        except Exception as e:
            raise RuntimeError(f'error opening genome file @ {genome_fp}') from e
    
        if len(ann_fps) != 2:
            raise ValueError(f'exactly two haplotypes supported')
        
        self.hindex = dict()
        self.hindex['h0'] = ann_fps[0]
        self.hindex['h1'] = ann_fps[1]
        self.gene_orders = dict()
        
        os.makedirs(out_dir, exist_ok = True)
        os.makedirs(temp_dir, exist_ok = True)
        self.out_dir = out_dir
        self.temp_dir = temp_dir
        self.min_pident = min_pident
        self.min_dlen = min_diag_len
        self.pad_len = pad_len
    
    def build_2d_matrix(
        self
    ):
        """
        args :
        returns :
        raises :
        """

        for hid, ann_fp in self.hindex.items():
            gene_order = utl.parse_protein_coding_genes(ann_fp)
            self.gene_orders[hid] = gene_order
        
        x = self.gene_orders["h0"]
        y = self.gene_orders["h1"]
        n = max(len(x), len(y))

        self.x_gorder = x
        self.y_gorder = y
        mat = np.zeros((n, n), dtype=int)

        for i in range(len(x)):
            g1, _, _, g1_str = x[i]
            for j in range(len(y)):
                g2, _, _, g2_str = y[j]
                if g1 == g2:
                    if g1_str == g2_str:
                        mat[i][j] = 1.0 # fwd match
                    else:
                        mat[i][j] = 2.0 # rev match
        return mat, n
    
    def extend_2d_matrix(
        self
    ):
        """
        """
        return

    







