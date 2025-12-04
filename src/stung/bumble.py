import pyfastx
import stung.utils as utl
import os
import numpy as np
from tqdm import tqdm
import subprocess

class Point:
    """
    A class representing a point in 2d space
    """
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
    """
    A class representing a block defined by two points, start and end.
    The start point is lower left corner; the end point is upper right corner.
    """
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
        return f'{str(self.start)}\t{str(self.end)}'

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
    write list of blocks into a text file
    """
    with open(out_fp, 'w') as fh:
        for b in blocks:
            fh.write(f'{str(b)}\n')

class Gene:
    def __init__(
        self,
        name : str,
        chr : str,
        start : int,
        end : int,
        strand : str
    ):
        """
        A class representing a gene.
        Retains a subset of properties that mjol offers for GFeature
        """
        self.name = name
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
    
    def __str__(
        self
    ):
        return f'{self.name}\t{self.chr}\t{self.start}\t{self.end}\t{self.strand}'

class Bumble:
    """
    A class representing the initial parameters for a stung analysis
    """
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
        initializes a Bumble instance
        args :
            self
            genome_fp (str) : genome FASTA file path
            ann_fps (list(str)) : list of file paths for gene annotation files (must be either gff/gtf formatted)
            temp_dir (str) : path to which intermediate files are saved
            out_dir (str) : output directory where global results (i.e., not specific to a breakpoint)
            min_pident (float; default = 90.0) : minimum percent identity to be considered a match
            min_diag_len (int; default = 10) : minimum diagonal length to be considered a colinear block
            pad_len (int; default = 10) : extra number of genes to analyze during colinear block refinement
        returns :
            None
        raises :
            None
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

        self.x_gene_order = None
        self.y_gene_order = None
        self.pident_mat = None
    
    def build_2d_matrix(self):
        """
        builds a 2d matrix based on HUGO gene symbol matches between
        protein coding genes in input annotations
        args :
            self
        returns :
            mat (np.ndarray) : square 2d matrix storing gene symbol matches
            n (int) : length of the matrix
        raises :
            None
        """

        for hid, ann_fp in self.hindex.items():
            gene_order, _ = utl.parse_protein_coding_genes(ann_fp)
            self.gene_orders[hid] = [
                Gene(
                    name = x[0],
                    chr = x[1],
                    start = x[2],
                    end = x[3],
                    strand = x[4]
                ) for x in gene_order
            ]
            self.gene_orders[hid].sort(key=lambda x: x.start)
        
        x = self.gene_orders["h0"]
        y = self.gene_orders["h1"]

        # write out gene_orders
        fn = os.path.join(self.out_dir, 'h0_gene_order.csv')
        with open(fn, 'w') as fh:
            for i, gene in enumerate(x):
                fh.write(f'{i},{gene.name}\n')
        
        fn = os.path.join(self.out_dir, 'h1_gene_order.csv')
        with open(fn, 'w') as fh:
            for i, gene in enumerate(y):
                fh.write(f'{i},{gene.name}\n')

        n = max(len(x), len(y))

        mat = np.zeros((n, n), dtype=int)

        for i in range(len(y)): # h1 on the y-axis
            g1 = y[i]
            for j in range(len(x)): # h0 on the x-axis
                g0 = x[j]
                if g1.name == g0.name:
                    if g1.strand == g0.strand:
                        mat[i][j] = 1.0 # same strand match
                    else:
                        mat[i][j] = 2.0 # opposite strand match
        self.x_gene_order = x
        self.y_gene_order = y

        return mat, n
    
    def extend_colinear_blocks(
        self,
        mat : np.ndarray,
        n : int,
        stung : tuple[Point, Point],
        bp_i : int
    ):
        """
        extends colinear blocks via pairwise sequence alignments that may add
        additional matches to the 2d match matrix
        args :
            mat (np.ndarray) : 2d match matrix
            n (int) : length of mat
            stung (tuple(Point, Point)) : a stung region defined by two points; left-closed, right-open (=lcro)
            bp_i (int) : breakpoint index
        returns :
            xlim (tuple(int, int)) : x axis limits (left and right); lcro
            ylim (tuple(int, int)) : y axis limits (top and bottom); lcro
            wkdir (str) : path to which intermediate files are saved
        raises :
            None
        """

        if self.pident_mat is None:
            self.pident_mat= np.zeros((n, n), dtype=float)

        start, end = stung

        xlim = (start.x - self.pad_len, end.x + self.pad_len)
        ylim = (start.y - self.pad_len, end.y + self.pad_len)

        wkdir = os.path.join(self.temp_dir, str(bp_i))
        cmd = f'mkdir -p {wkdir}'
        subprocess.call(cmd, shell=True)

        if not self.x_gene_order or not self.y_gene_order:
            raise ValueError()
        
        x = self.x_gene_order
        y = self.y_gene_order

        utl.plot_2d_matrix(
            mat = mat[ylim[0]:ylim[1], xlim[0]:xlim[1]],
            x_gene_order = None,
            y_gene_order = None,
            n = None,
            is_interactive = False,
            dot_s = 50, # bigger dots used
            save = True,
            out_fp = os.path.join(wkdir, f'bp_{bp_i}_pre.png'),
            x_offset = xlim[0],
            y_offset = ylim[0]
        )

        for i in tqdm(range(ylim[0], ylim[1], 1)): # row index
            g1 = y[i]

            s1 = get_genomic_seq(
                genome=self.genome,
                chr=g1.chr,
                start=g1.start,
                end=g1.end,
                strand=g1.strand
            )

            for j in range(xlim[0], xlim[1], 1): # iterates through col index

                # align two gene sequences if no match is found
                if mat[i][j] < 1:
                    g0 = x[j]

                    s0 = get_genomic_seq(
                        genome = self.genome,
                        chr=g0.chr,
                        start=g0.start,
                        end=g0.end,
                        strand=g0.strand
                    )

                    seq_fn = os.path.join(wkdir, 'to_aln.fa')
                    with open(seq_fn, 'w') as seq_fh:
                        seq_fh.write(f'>{g1.name}_{g1.chr}_{i}\n{s1.upper()}\n') 
                        seq_fh.write(f'>{g0.name}_{g0.chr}_{j}\n{s0.upper()}\n')

                    # run a*pa2 (messages suppressed)
                    aln_fn = os.path.join(wkdir, 'pa_aln.csv')
                    cmd = f'pa-bin --input {seq_fn} -o {aln_fn}'
                    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                
                    if os.path.getsize(aln_fn) == 0:
                        continue
                    
                    pident = utl.parse_pa_aln(aln_fn)
                    if pident > self.min_pident:
                        if g1.strand == g0.strand: # same strand match
                            mat[i][j] = 1.0
                        else: # opposite strand match
                            mat[i][j] = 2.0
                    self.pident_mat[i][j] = pident

                    cmd = f'rm {seq_fn} {aln_fn}'
                    subprocess.call(cmd, shell=True)
    
        utl.plot_2d_matrix(
            mat = mat[ylim[0]:ylim[1], xlim[0]:xlim[1]],
            x_gene_order = None,
            y_gene_order = None,
            n = None,
            is_interactive = False,
            dot_s = 50, # bigger dots used
            save = True,
            out_fp = os.path.join(wkdir, f'bp_{bp_i}_post.png'),
            x_offset = xlim[0],
            y_offset = ylim[0]
        )

        return xlim, ylim, wkdir

    def save_pident_mat(
        self,
        file_path : str,
        xlim : tuple[int, int],
        ylim : tuple[int, int]
    ):
        """
        saves pairwise percent identity matrix for a specific region
        args :
            file_path (str) : file path to which percent identity matrix is saved
            xlim (tuple(int, int)) : x axis limits (left, right); lcro
            ylim (tuple(int, int)) : y axis limits (top, bottom); lcro
        returns :
            None
        raises :
            None
        """
        with open(file_path, 'w') as fh:
            fh.write(f'i')
            for j in range(ylim[0], ylim[1]):
                fh.write(f'\t{j}')
            fh.write('\n')

            for i in range(xlim[0], xlim[1]):
                fh.write(f'{i}')
                for j in range(ylim[0], ylim[1]):
                    fh.write(f'\t{self.pident_mat[i][j]}')
                fh.write(f'\n')

def get_genomic_seq(
    genome : pyfastx.Fasta,
    chr : str,
    start : int,
    end : int,
    strand : str
) -> str:
    """
    extracts a genomic sequence from pyfastx.Fasta object
    args :
        chr (str) : chromosome
        start (int) : 0-based start coordinate
        end (int) : 0-based end coordinate
        strand (str) : strand information; either '+' or '-'
    returns :
        s (str) : extracted genomic sequence; upper-cased
    raises :
        ValueError
    """
    if strand != ['+', '-']:
        raise ValueError(f'unknown strand {strand}; must either be "+" or "-"')
    
    if strand == '+':
        s = genome[chr][start:end].seq.upper()
    else:
        s = genome[chr][start:end].antisense.upper()
    return s