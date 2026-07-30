#!/usr/bin/env python

from stung import bumble, stung, utils
import argparse
import os
import pandas as pd

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', type=str, help="csv", required=True)
    parser.add_argument('-g', '--genome', type=str, help="fasta", required=True)
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
    
        utils.save_2d_matrix(mat, os.path.join(wkdir, 'mmat.pre.tsv'))

        utils.plot_2d_matrix(
            mat = mat,
            pident_mat = None,
            x_gene_order = bmbl.x_gene_order,
            y_gene_order=bmbl.y_gene_order,
            n = n,
            is_interactive=False,
            save=True,
            out_file_path = os.path.join(wkdir, 'full.pre.png')
        )

        _, _, _ = stung.buzz(
            mat = mat,
            bmbl = bmbl,
            n = n,
            out_dir = wkdir,
            verbose = True
        )

        utils.save_2d_matrix(mat, os.path.join(wkdir, 'mmat.post.tsv'))

        utils.plot_2d_matrix(
            mat = mat,
            pident_mat = None,
            x_gene_order = bmbl.x_gene_order,
            y_gene_order=bmbl.y_gene_order,
            n = n,
            is_interactive=False,
            save=True,
            out_file_path = os.path.join(wkdir, 'full.post.png')
        )

        # NOTE: if not stung is detected in the first place, no pident_mat is constructed (i.e., pairwise alignment not executed)
        if bmbl.pident_mat is None:
            print(f"percent identity matrix not present; skipping the save() operation")
        else:
            utils.save_2d_matrix(bmbl.pident_mat, os.path.join(wkdir, 'pident.tsv'))

if __name__ == "__main__":
    main()