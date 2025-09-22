from mjol import gan
import os
import plotly.graph_objects as pgo
import numpy as np

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

            gene_order.append((gene_name, f.start, f.end, f.strand))
    return gene_order

def plot_2d_matrix(
    mat,
    gene_orders : dict,
    n : int,
    is_interactive : bool = True,
    save : bool = False
):
    """
    TODO
    """

    if len(gene_orders) == 0:
        raise ValueError()
    
    custom_data = np.empty((n, n), dtype=object)
    h0 = gene_orders['h0']
    h1 = gene_orders['h1']
    for i in range(n):
        h0_val = None if i >= len(h0) else h0[i][0]
        for j in range(n):
            h1_val = None if j >= len(h1) else h1[j][0]
            custom_data[i][j] = (h0_val, h1_val)
    
    if is_interactive:
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
        # TODO: finish implementation
        raise NotImplementedError