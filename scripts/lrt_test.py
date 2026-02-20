#!/usr/bin/env python

"""
author :  HJ Ji
date : 02/20/26
contact : hji20@jh.edu
"""

import argparse
from Bio import Phylo
from io import StringIO
from scipy import stats
import sys

def calculate_lrt(
        tree_str : str, 
        lnL_m0 : float,
        lnL_m1 : float,
        cutoff : float
    ):
    tree = Phylo.read(StringIO(tree_str), "newick")
    num_taxa = tree.count_terminals()
    num_branches = (2 * num_taxa) - 3
    df = num_branches - 1 # number of params ~ ratios
    lrt_stat = 2 * (lnL_m1 - lnL_m0)
    pval = stats.chi2.sf(lrt_stat, df)
    
    print("-" * 30, file=sys.stderr)
    print(f"Tree Analysis Results", file=sys.stderr)
    print("-" * 30, file=sys.stderr)
    print(f"Total Branches (Parameters in M1): {num_branches}", file=sys.stderr)
    print(f"Degrees of Freedom (df):           {df}", file=sys.stderr)
    print(f"LRT Statistic (2ΔlnL):            {lrt_stat:.4f}", file=sys.stderr)
    print(f"P-value:                          {pval:.4e}", file=sys.stderr)
    print("-" * 30, file=sys.stderr)

    if pval < cutoff:
        print("SIGNIFICANT: Reject m0 in favor of m1 (Free-ratio).", file=sys.stderr)
    else:
        print("NOT SIGNIFICANT: Fail to reject m0. Use the simpler model.", file=sys.stderr)

    print(pval)

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-tree", type=str, help="", required=True)
    parser.add_argument("-m0", type=float, help="", required=True)
    parser.add_argument("-m1", type=float, help="", required=True)
    parser.add_argument("-c", "--cut-off", type=float, help="", required=False, default=0.05)
    
    args = parser.parse_args()
    calculate_lrt(args.tree, args.m0, args.m1, args.cut_off)

if __name__ == "__main__":
    main()