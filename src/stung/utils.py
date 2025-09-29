from mjol import gan
import os
import plotly.graph_objects as pgo
import matplotlib.pyplot as plt
import re
import numpy as np

EXPECTED_PLT_DIM=100

def is_protein_coding(
    gene : gan.GFeature
) -> bool:
    """
    """
    for tx in gene.children:
        for child in tx.children:
            if child.feature_type == "CDS":
                return True
    return False

def parse_protein_coding_genes(
    fp : str
):
    """
    """
    fmt = os.path.basename(fp).split('.')[-1].lower()
    if fmt not in ['gff', 'gtf']:
        raise ValueError()
    
    try:
        gan_db = gan.GAn(
            file_name = fp,
            file_fmt = fmt
        )
        gan_db.build_db()
    except Exception as e:
        raise RuntimeError()
    
    gene_order = []
    for uid in gan_db.features:
        f = gan_db.get_feature(
            uid = uid
        )
        if f.feature_type == "gene" and is_protein_coding(f):
            gene_name = f.attributes['gene_name'] if 'gene_name' in f.attributes else None
            
            if not gene_name:
                print("warning"); continue

            gene_order.append((gene_name, f.chr, f.start, f.end, f.strand))
    return gene_order

def plot_2d_matrix(
    mat,
    x_gene_order : list,
    y_gene_order : list,
    n : int,
    is_interactive : bool = True,
    fig_size : tuple = (10, 10),
    dpi : int = 100,
    dot_s : int = 1,
    save : bool = False,
    out_fp : str = None,
    x_offset : int = 0,
    y_offset : int = 0
):
    """
    TODO
    """
    
    if save and not out_fp:
        raise ValueError()
    
    if save and is_interactive:
        raise ValueError()
    
    if (x_offset > 0 or y_offset > 0) and is_interactive:
        raise ValueError()
    
    if is_interactive:
        if n != len(mat) or n != len(mat[0]):
            raise ValueError()
        
        custom_data = np.empty((n, n), dtype=object)
        for i in range(n):
            for j in range(n):
                if mat[i][j] == 0:
                    custom_data[i][j] = ("None", "None")
                else:
                    y_val = None if i >= len(y_gene_order) else y_gene_order[i].name
                    x_val = None if j >= len(x_gene_order) else x_gene_order[j].name
                    custom_data[i][j] = (x_val, y_val)
        custom_colorscale = [
            [0.0, "white"],
            [0.333, "white"],
            [0.334, "blue"],
            [0.666, "blue"],
            [0.667, "red"],
            [1.0, "red"]
        ]
        fig = pgo.Figure(data=pgo.Heatmap(
            z=mat,
            colorscale=custom_colorscale,
            customdata=custom_data,
            showscale=False,
            hovertemplate="x: %{customdata[0]}<br>y: %{customdata[1]}<extra></extra>",
        ))

        fig.update_layout(
            dragmode='zoom',
            xaxis=dict(
                scaleanchor='y',
                range=[0, n]
            ),
            yaxis = dict(
                range=[0, n]
            )
        )

        fig.show()
    else:
        y_1, x_1 = np.where(mat == 1) # fwd matches
        y_2, x_2 = np.where(mat == 2) # rev matches

        fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)

        ax.set_ylabel('h1')
        ax.set_xlabel('h0')

        ax.scatter(x_1, y_1, color='blue', marker='s', s=dot_s)
        ax.scatter(x_2, y_2, color='red', marker='s', s=dot_s)

        ax.set_xlim(0, len(mat[0]))
        ax.set_ylim(0, len(mat))

        if x_offset > 0:
            if len(mat[0]) > EXPECTED_PLT_DIM:
                print(f'warning : number of cols ({len(mat[0])}) may be too large for vis')
            x_ticks = np.arange(len(mat[0]))
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_ticks + x_offset, rotation=90)
        
        if y_offset > 0:
            if len(mat) > EXPECTED_PLT_DIM:
                print(f'warning : number of rows {len(mat)}) may be too large for vis')
            y_ticks = np.arange(len(mat))
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_ticks + y_offset, rotation=90)

        if save:
            fig.savefig(out_fp, bbox_inches='tight')
        else:
            plt.show()

        plt.close(fig)


def parse_cigar(
    s : str
) -> list[tuple[str, int]]:
    """
    TODO
    """
    cigar_ops = []
    i = 0
    while i < len(s):
        res = re.match(r'\d+', s[i:])
        if res:
            match_s = res.group()
            op_len = int(match_s)
            i += len(match_s)
            op_char = s[i]
            cigar_ops.append((op_len, op_char))
            i += 1
        else:
            cigar_ops.append((1, s[i]))
            i += 1
    
    return cigar_ops

def parse_pa_aln(
    fp : str
):
    """
    TODO
    """
    pident = None
    with open(fp, 'r') as fh:
        for ln in fh:
            parts = ln.strip().split(",")
            cigar = parts[1]
            cigar_ops = parse_cigar(cigar)
            tot_l = sum([x[0] for x in cigar_ops])
            eq_l = sum([x[0] for x in cigar_ops if x[1] == '='])
            pident = eq_l / tot_l * 100 # TODO: discuss this calculation
    return pident