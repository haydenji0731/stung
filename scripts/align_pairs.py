#!/usr/bin/env python

"""
author : HJ Ji
date : 02/02/26
contact : hji20@jh.edu
"""

import argparse
import os
from stung import utils
import pyfastx
import subprocess
import numpy as np

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-h0', '--h0-file', type=str, help="", required=True)
    parser.add_argument('-h1', '--h1-file', type=str, help="", required=True)
    parser.add_argument('-q', '--query-gene', type=str, help="", required=True)
    parser.add_argument('-g', '--genome', type=str, help="", required=True)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=False, default='out')
    parser.add_argument('-x', '--x-coords', type=int, nargs=2, help="fully closed gene order interval", required=True)
    parser.add_argument('-y', '--y-coords', type=int, nargs=2, help="fully closed gene order interval", required=True)
    args = parser.parse_args()
    return args

def check_args(args) -> str:
    if not os.path.exists(args.h0_file) or \
        not os.path.exists(args.h1_file) or \
        not os.path.exists(args.genome):
        raise FileNotFoundError(f'input file(s) dne')
    os.makedirs(args.out_dir, exist_ok=True)
    temp_dir = os.path.join(args.out_dir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir

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
    query_gene : str,
    x_chrom : str,
    x_gene_order : list,
    x_coords : tuple[int, int],
    y_chrom : str,
    y_gene_order : list,
    y_coords : tuple[int, int],
    out_dir : str
):
    xst, xen = x_coords
    x_chr = genome[x_chrom].seq
    xi_set = set()
    yi_set = set()
    for xi in range(xst, xen + 1, 1):
        x_gene = x_gene_order[xi]
        if query_gene not in x_gene[0]:
            continue
        xi_set.add(xi)
        file_name = os.path.join(out_dir, f'x_{xi}.fa')
        gst = x_gene[2] - 1 # adjust to 0-based, half-closed, half-open interval
        gen = x_gene[3]
        s = x_chr[gst:gen]
        with open(file_name, 'w') as file_handle:
            file_handle.write(f'>x_{xi}\n{s}\n')

    yst, yen = y_coords
    y_chr = genome[y_chrom].seq
    for yi in range(yst, yen + 1, 1):
        y_gene = y_gene_order[yi]
        if query_gene not in y_gene[0]:
            continue
        yi_set.add(yi)
        file_name = os.path.join(out_dir, f'y_{yi}.fa')
        gst = y_gene[2] - 1 # adjust to 0-based, half-closed, half-open interval
        gen = y_gene[3]
        s = y_chr[gst:gen]
        with open(file_name, 'w') as file_handle:
            file_handle.write(f'>y_{yi}\n{s}\n')
    return xi_set, yi_set

def align(
    temp_dir : str,
    x_coords : tuple[int, int],
    y_coords : tuple[int, int],
    xi_set : set[int],
    yi_set : set[int]
) -> np.ndarray:
    xst, xen = x_coords
    yst, yen = y_coords
    width = xen - xst + 1
    height = yen - yst + 1

    res = np.zeros((height, width), dtype=np.float64)

    for xi in range(xst, xen + 1, 1):
        if xi not in xi_set:
            continue
        x_file_name = os.path.join(temp_dir, f'x_{xi}.fa')
        for yi in range(yst, yen + 1, 1):
            if yi not in yi_set:
                continue
            y_file_name = os.path.join(temp_dir, f'y_{yi}.fa')
            seq_fn = os.path.join(temp_dir, 'to_aln.fa')
            cmd = f'cat {x_file_name} {y_file_name} > {seq_fn}'
            subprocess.call(cmd, shell=True)

            # run a*pa2 (messages suppressed)
            aln_fn = os.path.join(temp_dir, 'pa_aln.csv')
            cmd = f'pa-bin --input {seq_fn} -o {aln_fn}'
            subprocess.call(
                cmd, shell=True, stdout=subprocess.DEVNULL, 
                stderr=subprocess.DEVNULL
            )

            if os.path.getsize(aln_fn) == 0:
                continue
        
            pident = utils.parse_pa_aln(aln_fn)
            row_i = yi - yst
            col_i = xi - xst
            res[row_i, col_i] = pident

            cmd = f'rm {seq_fn} {aln_fn}'
            subprocess.call(cmd, shell=True)
    return res

def main():
    args = parse()
    temp_dir = check_args(args)

    print(f'INFO - parsing h0 and h1 annotations')
    h0_gene_order, _ = utils.parse_protein_coding_genes(args.h0_file)
    h1_gene_order, _ = utils.parse_protein_coding_genes(args.h1_file)

    x_coords = args.x_coords
    y_coords = args.y_coords

    check_coords(
        x_coords = x_coords,
        y_coords = y_coords,
        h0_len = len(h0_gene_order),
        h1_len = len(h1_gene_order)
    )
    
    print(f'INFO - extracting gene sequences (introns included)')
    genome = pyfastx.Fasta(args.genome)
    xi_set, yi_set = extract(
        genome = genome,
        query_gene = args.query_gene,
        x_chrom = h0_gene_order[0][1],
        x_gene_order = h0_gene_order,
        x_coords = x_coords,
        y_chrom = h1_gene_order[0][1],
        y_gene_order = h1_gene_order,
        y_coords = y_coords,
        out_dir = temp_dir
    )

    print(f'INFO - running pairwise alignments')
    res = align(
        temp_dir = temp_dir,
        x_coords = x_coords,
        y_coords = y_coords,
        xi_set = xi_set,
        yi_set = yi_set
    )

    xst, xen = x_coords
    yst, yen = y_coords
    out_file_path = os.path.join(args.out_dir, f'x_{xst}:{xen}_y_{yst}:{yen}.pident.tsv')
    utils.save_2d_matrix(res, out_file_path)
    print(f'INFO - Done!')

if __name__ == "__main__":
    main()



    