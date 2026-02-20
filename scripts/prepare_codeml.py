#!/usr/bin/env python

"""
author :  HJ Ji
date : 02/18/26
contact : hji20@jh.edu
"""

import os
import pandas as pd
import pyfastx
import argparse
from collections import defaultdict
import re

STOP_CODONS = ['taa', 'tag', 'tga']

# positive int type def'n
def pint(val):
    ival = int(val)
    if ival <= 0:
        raise argparse.ArgumentTypeError(f"{val} is not a positive integer")
    return ival

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-p', '--protein-file', type=str, help="", required=True)
    parser.add_argument('-x', '--cds-file', type=str, help="", required=True)
    parser.add_argument('-g', '--genes-file', type=str, help="", required=True)
    parser.add_argument('-l', '--list-file',  type=str, help="", required=True)
    parser.add_argument('-c', '--chrom', type=str, help="", required=True)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=False, default='.')
    args = parser.parse_args()
    return args

def check_args(args):
    if not os.path.exists(args.protein_file) or \
        not os.path.exists(args.cds_file) or \
        not os.path.exists(args.genes_file) or \
        not os.path.exists(args.list_file):
        raise FileNotFoundError(f'input file(s) dne')
    # TODO: add a check for chrom?
    os.makedirs(args.out_dir, exist_ok=True)

def load_genes_df(file_path : str):
    """
    loads ordered list of genes from a CSV file
    """
    df = pd.read_csv(file_path, header=None)
    df.columns = [
        'index', 'gene_name', 'gene_id', 'start', 'end', 'strand'
    ]
    return df

def load_list_of_indices(file_path : str) -> list[int]:
    indices = []
    with open(file_path, 'r') as file:
        for ln in file:
            indices.append(int(ln.strip()))
    return indices

def select_tx(
    gene_id : str,
    protein_id_lookup : dict,
    proteins : pyfastx.Fasta
) -> str:
    """
    selects the transcript encoding longest protein for each gene
    """
    if gene_id not in protein_id_lookup:
        raise KeyError(f'{gene_id} not found in gene_lookup')
    
    plens = []

    for tx in protein_id_lookup[gene_id]:
        if tx not in proteins:
            raise KeyError(f'{tx} not found in list of proteins')
        pseq = proteins[tx].seq
        if '.' in pseq: # remove stop codons
            pseq_mod = pseq.replace('.', '')
            plens.append((tx, len(pseq_mod)))
        else:
            plens.append((tx, len(pseq)))
    
    return max(plens, key=lambda x : x[1])

def build_protein_lookup(proteins):
    protein_id_lookup = defaultdict(list)
    for x in proteins:
        name = x.name
        protein_id_lookup[name.split('.')[0]].append(name)
    return protein_id_lookup

def write(
    out_dir : str,
    txes : dict[str],
    gene_name_lookup : str,
    proteins : pyfastx.Fasta,
    cdses : pyfastx.Fasta
):
    out_file_path1 = os.path.join(out_dir, 'proteins.fa')
    out_file_path2 = os.path.join(out_dir, 'translated_nucls.fa')

    ctr = defaultdict(int)

    outfile1 = open(out_file_path1, 'w')
    outfile2 = open(out_file_path2, 'w')

    for gene, tx in txes.items():
        gene_name = gene_name_lookup[gene]
        gi = ctr[gene_name]

        pseq = proteins[tx].seq
        xseq = cdses[tx].seq
        res = [m.start() for m in re.finditer(r'\.', pseq)]

        pseq = pseq.replace('.', '')
        outfile1.write(f'>{gene_name}_{gi}\n{pseq}\n')

        if res:
            xseq_mod = ""
            prev_pos = None
            for i in range(0, len(res), 1):
                pos = res[i] * 3
                if not prev_pos:
                    xseq_mod += xseq[:pos]
                else:
                    xseq_mod += xseq[prev_pos:pos]
                prev_pos = pos + 3
            xseq_mod += xseq[prev_pos:]
            xseq = xseq_mod
        
        last_codon = xseq[-3:].lower()
        if last_codon in STOP_CODONS:
            xseq = xseq[:-3]
        
        outfile2.write(f'>{gene_name}_{gi}\n{xseq}\n')

        ctr[gene_name] += 1

    outfile1.close()
    outfile2.close()

def main():
    args = parse()
    check_args(args)

    try:
        cdses = pyfastx.Fasta(args.cds_file)
        proteins = pyfastx.Fasta(args.protein_file)
    except Exception as e:
        raise RuntimeError(f'error while open cds and/or protein file : {e}')

    protein_id_lookup = build_protein_lookup(proteins)
    genes = load_genes_df(args.genes_file)
    gene_name_lookup = dict(zip(genes['gene_id'], genes['gene_name']))
    indices = load_list_of_indices(args.list_file)

    try:
        qry_genes = genes.iloc[indices]["gene_id"].to_list()
    except Exception as e:
        raise RuntimeError(f'error subsetting from genes_df : {e}')
    
    selected_txes = dict()
    for gene in qry_genes:
        res = select_tx(gene, protein_id_lookup, proteins)
        selected_txes[gene] = res[0]
    
    write(args.out_dir, selected_txes, gene_name_lookup, proteins, cdses)

if __name__ == "__main__":
    main()