#!/usr/bin/env python

"""
author : HJ Ji
date : 09/29/25
contact : hji20@jh.edu
"""

from stung import bumble, stung, utils
import argparse
import os
import pandas as pd

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', type=str, help="", required=True)
    parser.add_argument('-g', '--genome', type=str, help="", required=True)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=False, default='out')
    args = parser.parse_args()
    return args

def load_input(
    file_path : str
):
    df = pd.read_csv(file_path, header=None, comment='#')
    df.columns = [
        "prefix",
        "h0_filepath",
        "h1_filepath"
    ]
    return df
    
def main():
    args = parse()

    if not os.path.exists(args.input) or \
        not os.path.exists(args.genome):
        raise FileNotFoundError()

    os.makedirs(args.out_dir, exist_ok=True)

    in_df = load_input(args.input)

    for _, row in in_df.iterrows():
        prefix = row['prefix']
        h0_fp = row['h0_filepath']
        h1_fp = row['h1_filepath']
        print(f"currently processing...\nh0 : {h0_fp}\nh1 : {h1_fp}")
        
        if not os.path.exists(h0_fp) or \
            not os.path.exists(h1_fp):
            raise FileNotFoundError

        wkdir = os.path.join(args.out_dir, prefix)
        temp_dir = os.path.join(wkdir, 'temp')
        
        bmbl = bumble.Bumble(
            genome_fp = args.genome,
            ann_fps = [
                h0_fp,
                h1_fp
            ],
            temp_dir = temp_dir,
            out_dir = wkdir
        )
        mat, n = bmbl.build_2d_matrix()
        utils.save_2d_matrix(mat, os.path.join(wkdir, 'mat.pre.tsv'))

        utils.plot_2d_matrix(
            mat = mat,
            x_gene_order = bmbl.x_gene_order,
            y_gene_order=bmbl.y_gene_order,
            n = n,
            is_interactive=False,
            save=True,
            out_fp = os.path.join(wkdir, 'full.pre.png')
        )

        # TODO : is this optimal?
        _, _, _ = stung.buzz(
            mat = mat,
            bmbl = bmbl,
            n = n,
            out_dir = wkdir,
            verbose = True
        )
        utils.save_2d_matrix(mat, os.path.join(wkdir, 'mat.post.tsv'))

        utils.plot_2d_matrix(
            mat = mat,
            x_gene_order = bmbl.x_gene_order,
            y_gene_order=bmbl.y_gene_order,
            n = n,
            is_interactive=False,
            save=True,
            out_fp = os.path.join(wkdir, 'full.post.png')
        )

if __name__ == "__main__":
    main()
