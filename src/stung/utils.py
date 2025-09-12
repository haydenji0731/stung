from mjol import gan
import os
import plotly.graph_objects as pgo

def is_protein_coding(
    gene : gan.gFeature
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
    gene_orders,
    is_interactive : bool = False,
    save : bool = False
):
    """
    """
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
            showscale=False,  # hide colorbar
            hovertemplate="x: %{x}<br>y: %{y}<br>value: %{z}<extra></extra>",
        ))

        fig.update_layout(
            dragmode='zoom',
            xaxis=dict(scaleanchor='y'),
        )

        fig.show()
    else:
        raise NotImplementedError