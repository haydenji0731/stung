#!/usr/bin/env python

"""
author :  HJ Ji
date : 10/02/25
contact : hji20@jh.edu
"""

# dotplotter available @ https://github.com/drboothtj/dotplotter/tree/main

import argparse
import os
from stung import utils
import pyfastx
import subprocess

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', '--blastn-path', type=str, help="", required=False, default='blastn')
    parser.add_argument('-h0', '--h0-file', type=str, help="", required=True)
    parser.add_argument('-h1', '--h1-file', type=str, help="", required=True)
    parser.add_argument('-g', '--genome', type=str, help="", required=True)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=False, default='out')
    #parser.add_argument('-f', '--plt-format', choices=["png", "svg"], required=False, default='svg')
    parser.add_argument('-x', '--x-coords', type=int, nargs=2, help="fully closed", required=True)
    parser.add_argument('-y', '--y-coords', type=int, nargs=2, help="fully closed", required=True)
    args = parser.parse_args()
    return args

def check_args(args):
    if not os.path.exists(args.h0_file) or \
        not os.path.exists(args.h1_file) or \
        not os.path.exists(args.genome):
        raise FileNotFoundError(f'input file(s) dne')
    os.makedirs(args.out_dir, exist_ok=True)

def check_coords(
    x_coords: tuple[int, int],
    y_coords: tuple[int, int],
    h0_len: int,
    h1_len: int
):
    x0, x1 = x_coords
    y0, y1 = y_coords

    if any(coord < 0 for coord in (x0, x1, y0, y1)):
        raise ValueError("coordinates must be non-negative")

    if not (x0 < x1 < h0_len):
        raise ValueError("invalid x range")

    if not (y0 < y1 < h1_len):
        raise ValueError("invalid y range")
    
def extract(
    genome : pyfastx.Fasta,
    x_chrom : str,
    x_grange : tuple[int, int],
    y_chrom : str,
    y_grange : tuple[int, int],
    out_dir : str
):
    """
    extracts genomic sequences to align
    """
    if x_chrom not in genome or y_chrom not in genome:
        raise KeyError()
    
    x_gstart, x_gend = x_grange
    y_gstart, y_gend = y_grange
    x_gseq = genome[x_chrom].seq[x_gstart:x_gend]
    y_gseq = genome[y_chrom].seq[y_gstart:y_gend]
    
    x_out_fn = os.path.join(out_dir, 'x_gseq.fa')
    with open(x_out_fn, 'w') as fh:
        fh.write(f'>{x_chrom}_{x_gstart}:{x_gend}\n{x_gseq}')

    # fai index helps compare the sequence lengths
    cmd = f'samtools faidx {x_out_fn}'
    subprocess.call(cmd, shell=True)
    
    y_out_fn = os.path.join(out_dir, 'y_gseq.fa')
    with open(y_out_fn, 'w') as fh:
        fh.write(f'>{y_chrom}_{y_gstart}:{y_gend}\n{y_gseq}')
    cmd = f'samtools faidx {y_out_fn}'
    subprocess.call(cmd, shell=True)

    return x_out_fn, y_out_fn

def align(
    x_fn : str,
    y_fn : str,
    blastn_path : str,
    aln_fn : str
):
    """
    run blastn alignments
    """
    if not os.path.exists(x_fn) or not os.path.exists(y_fn):
        raise FileNotFoundError()
    
    cmd = f'{blastn_path} -query {x_fn} -subject {y_fn} -outfmt 6 > {aln_fn}'
    subprocess.call(cmd, shell=True)
    
def main():
    args = parse()
    check_args(args)

    print(f'INFO - parsing h0 and h1 annotations')
    h0_gene_order, _ = utils.parse_protein_coding_genes(args.h0_file)
    h1_gene_order, _ = utils.parse_protein_coding_genes(args.h1_file)

    check_coords(
        x_coords = args.x_coords,
        y_coords = args.y_coords,
        h0_len = len(h0_gene_order),
        h1_len = len(h1_gene_order)
    )

    genome = pyfastx.Fasta(args.genome)
    
    x_start, x_end = args.x_coords
    x_gstart = h0_gene_order[x_start][3] - 1 # adjust to 0-based, half-closed, half-open interval
    x_gend = h0_gene_order[x_end][4]
    y_start, y_end = args.y_coords
    y_gstart = h1_gene_order[y_start][3] - 1 # adjust to 0-based, half-closed, half-open interval
    y_gend = h1_gene_order[y_end][4]

    print(f'query   [{x_start}:{x_end}]  {x_gstart:,} ~ {x_gend:,}')
    print(f'subject [{y_start}:{y_end}]  {y_gstart:,} ~ {y_gend:,}')

    print(f'INFO - extracting genomic sequences')
    x_fn, y_fn = extract(
        genome = genome,
        x_chrom = h0_gene_order[0][2],
        x_grange = (x_gstart, x_gend),
        y_chrom = h1_gene_order[0][2],
        y_grange = (y_gstart, y_gend),
        out_dir = args.out_dir
    )

    print(f'INFO - running blastn')
    aln_fn = os.path.join(args.out_dir, f'x_{x_start}:{x_end}_y_{y_start}:{y_end}.blastn.tsv')
    align(
        x_fn = x_fn,
        y_fn = y_fn,
        blastn_path = args.blastn_path,
        aln_fn = aln_fn
    )

if __name__ == "__main__":
    main()