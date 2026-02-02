from mjol import gan
import os
import pandas as pd
import plotly.graph_objects as pgo
import matplotlib.pyplot as plt
import re
import numpy as np

EXPECTED_PLT_DIM=100

def is_protein_coding(
    gene : gan.GFeature
) -> bool:
    """
    checks if a GFeature instance is a protein coding gene
    based on whether or not the gene has a CDS record associated with it
    args :
        gene (gan.GFeature) : input gene as GFeature instance
    returns :
        is_protein_coding : True if the input gene is protein coding
    raises :
        None
    """
    is_protein_coding = False
    for tx in gene.children:
        for child in tx.children:
            if child.feature_type == "CDS":
                is_protein_coding = True
    return is_protein_coding

def parse_protein_coding_genes(
    file_path : str
):
    """
    parses a list of protein coding genes from an annotation file
    args :
        file_path (str) : path to an input gene annotation file
    returns :
        gene_order (list) : list of protein coding genes as tuples
        gan_db (gan.GAn) : gan database instance constructed from the annotation
    raises :
        ValueError
    """
    fmt = os.path.basename(file_path).split('.')[-1].lower()
    if fmt not in ['gff', 'gtf']:
        raise ValueError(f'invalid annotation file format {fmt}')
    
    try:
        gan_db = gan.GAn(
            file_name = file_path,
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
                print(f"warning : no gene name detected for feature {f}"); continue

            gene_order.append((gene_name, f.chr, f.start, f.end, f.strand, f))
    return gene_order, gan_db

# plotting helper functions

def do_labels_overlap(ax, axis='x') -> bool:
    fig = ax.figure
    fig.canvas.draw()
    
    if axis == 'x':
        labels = ax.get_xticklabels()
    else:
        labels = ax.get_yticklabels()

    bboxes = [label.get_window_extent() for label in labels if label.get_text()]

    for i, bbox1 in enumerate(bboxes):
        for bbox2 in bboxes[i+1:]:
            if bbox1.overlaps(bbox2):
                return True
    return False

def auto_resize_until_no_overlap(fig, ax, axis='x', max_iter=5, scale=1.2):
    for _ in range(max_iter):
        if not do_labels_overlap(ax, axis):
            return fig

        w, h = fig.get_size_inches()
        if axis == 'x':
            fig.set_size_inches((w * scale, h))
        else:
            fig.set_size_inches((w, h * scale))
        fig.canvas.draw()
    return fig

def load_2d_matrix(
    file_path : str
) -> np.ndarray:
    """
    TODO
    args :
    returns :
    raises :
        FileNotFoundError
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'file @ {file_path} dne')
    
    mat = pd.read_csv(file_path, sep=r"\s+", index_col=0)
    mat = mat.to_numpy()
    return mat

def plot_2d_matrix(
    mat : np.ndarray,
    pident_mat : np.ndarray = None,
    x_gene_order : list[str] = None,
    y_gene_order : list[str] = None,
    n : int = 0,
    is_interactive : bool = True,
    fig_size : tuple = (10, 10), # starting default plot size
    dpi : int = 100,
    dot_s : int = 1,
    save : bool = False,
    out_file_path : str = None,
    x_offset : int = 0,
    y_offset : int = 0
):
    """
    plots 2d matrix as a static or interactive plot
    args :
        mat (np.ndarray) : match matrix
        pident_mat (nd.ndarray) : pairwise percent identity matrix
        x_gene_order (list(str)) : ordered list of x-axis genes
        y_gene_order (list(str)) : ordered list of y-axis genes
        n (int) : match matrix dimension
        is_interactive (bool; default=True) : whether to enable interactive plotting
        fig_size (tuple(int, int); default=(10, 10)) : figure size for static plotting
        dpi (int; default=100) : dpi for static plotting
        dot_s (int; default=1) : dot size for static plotting
        save (bool; default=False) : whether to save the static plot
        out_file_path (str) : output file path where a plot is saved
        x_offset (int; default=0) : offset used to adjust the x-axis indices
        y_offset (int; default=0) : offset used to adjust the y-axis indices
    returns :
        None
    raises :
        ValueError
        RuntimeError
    """
    if is_interactive:

        # check args
        if save:
            raise ValueError('save operation not supported for interactive plots')
        
        if (x_offset > 0 or y_offset > 0):
            raise ValueError(f'x and/or y offsets not supported for interactive plots')
        
        if not x_gene_order or not y_gene_order:
            raise ValueError(f'please provide x and y-axis gene lists')
            
        if n != len(mat):
            raise ValueError(f'invalid dimension {n}; please provide the correct dimension for interactive plotting')
        
        custom_data = [[None for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if mat[i][j] == 0:
                    custom_data[i][j] = ("None", "None", "None")
                else:
                    y_val = None if i >= len(y_gene_order) else y_gene_order[i].name
                    x_val = None if j >= len(x_gene_order) else x_gene_order[j].name
                    pident = None if (pident_mat is None or pident_mat[i][j] == 0.0) else f"{pident_mat[i][j]:.2f}"
                    custom_data[i][j] = [x_val, y_val, str(pident)]
        custom_colorscale = [
            [0.0, "white"],  # 0 -> white
            [0.5, "blue"],   # 1 -> blue
            [1.0, "red"]     # 2 -> red
        ]
        fig = pgo.Figure(data=pgo.Heatmap(
            z=mat,
            colorscale=custom_colorscale,
            customdata=custom_data,
            showscale=False,
            zmin=0,
            zmax=2,
            hovertemplate=(
                "x: %{customdata[0]}<br>"
                "y: %{customdata[1]}<br>"
                "percent identity: %{customdata[2]}<br>"
                "<extra></extra>"
            )
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
            ax.set_yticklabels(y_ticks + y_offset)

        # automatic resizing
        fig = auto_resize_until_no_overlap(fig, ax, axis='x')
        fig = auto_resize_until_no_overlap(fig, ax, axis='y')

        if save:
            try:
                fig.savefig(out_file_path, bbox_inches='tight')
            except Exception as e:
                raise RuntimeError(f'error while saving the plot : {e}') from e
        else:
            plt.show()

        plt.close(fig)

class Region:
    """
    A class representing a subregion of a 2d matrix
    """
    def __init__(
        self,
        chr : str,
        mat : np.ndarray,
        x_start : int,
        x_end : int,
        y_start : int,
        y_end : int
    ):
        self.chr = chr
        self.xst = x_start
        self.xen = x_end
        self.yst = y_start
        self.yen = y_end
        self.data = mat[y_start:y_end, x_start:x_end].copy()
    
    def __str__(
        self
    ):
        return f"x:{self.xst}-{self.xen}\ny:{self.yst}-{self.yen}"

MAJOR_TICK_SIZE = 6
MINOR_TICK_SIZE = 3

def plot_2d_matrices(
    regions : list[Region],
    dot_sizes : list[int],
    out_file_path : str,
    fig_size : tuple = (3.42, 3.42),
    file_type : str = "svg",
    axis_title_size : int = 8,
    tick_label_size : int = 6,
    xlabel_vert : float = 0.01,
    ylabel_vert : float = 0.55,
):
    """
    TODO
    args :
    returns :
    raises :
        ValueError
    """
    if len(regions) % 2 != 0:
        raise ValueError(f"number of regions must be divisible by two; got {len(regions)}")

    if len(regions) != len(dot_sizes):
        raise ValueError()
    
    n_regions = len(regions)
    n_cols = 2
    n_rows = int(np.ceil(n_regions / n_cols))

    fig, axes = plt.subplots(
        nrows = n_rows,
        ncols = n_cols,
        figsize = (fig_size[0] * n_cols, fig_size[1] * n_rows)
    )
    
    axes = np.atleast_2d(axes)
    fig.patch.set_alpha(0)

    for i, region in enumerate(regions):
        dot_s = dot_sizes[i]
        x_offset = region.xst
        y_offset = region.yst

        y_1, x_1 = np.where(region.data == 1) # fwd matches
        y_2, x_2 = np.where(region.data == 2) # rev matches

        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        ax.patch.set_alpha(0)

        # draw dots
        ax.scatter(x_1, y_1, color="blue", marker="s", s=dot_s)
        ax.scatter(x_2, y_2, color="red", marker="s", s=dot_s)
        ax.set_xlim(0, region.data.shape[1])
        ax.set_ylim(0, region.data.shape[0])

        ax.tick_params(axis="both", labelsize=tick_label_size)

        if x_offset > 0:
            x_ticks = np.arange(region.data.shape[1] + 1)
            ax.set_xticks(x_ticks)

            # add tick label each increment of 5
            x_labels = [
                str(t + x_offset) if (t + x_offset) % 5 == 0 else ""
                for t in x_ticks
            ]

            ax.set_xticklabels(x_labels, rotation=90)

            for tick in ax.xaxis.get_ticklines():
                tick_val = tick.get_xdata()[0]
                if (tick_val + x_offset) % 5 == 0:
                    tick.set_markersize(MAJOR_TICK_SIZE)
                else:
                    tick.set_markersize(MINOR_TICK_SIZE)
            
        if y_offset > 0:
            y_ticks = np.arange(region.data.shape[0] + 1)
            ax.set_yticks(y_ticks)

            # add tick label each increment of 5
            y_labels = [
                str(t + y_offset) if (t + y_offset) % 5 == 0 else ""
                for t in y_ticks
            ]

            ax.set_yticklabels(y_labels, rotation=90)
            for tick in ax.yaxis.get_ticklines():
                tick_val = tick.get_ydata()[0]
                if (tick_val + y_offset) % 5 == 0:
                    tick.set_markersize(MAJOR_TICK_SIZE)
                else:
                    tick.set_markersize(MINOR_TICK_SIZE)
    
    # add superlabels to x and y axis
    fig.supxlabel("Maternal", fontsize=axis_title_size, y=xlabel_vert)
    fig.supylabel("Paternal", fontsize=axis_title_size, y=ylabel_vert)
    
    fig.tight_layout()
    fig.savefig(out_file_path, format=file_type, bbox_inches='tight')

def save_2d_matrix(
        mat : np.ndarray,
        file_path : str
):
    """
    saves a 2d matrix as a TSV file
    args :
        mat (np.ndarray) : input 2d matrix
        file_path (str) : where the matrix is saved to
    returns :
        None
    raises :
        None
    """
    n_cols = len(mat[0])
    n_rows = len(mat)
    with open(file_path, 'w') as fh:
        fh.write(f'i')
        for j in range(n_cols):
            fh.write(f'\t{j}')
        fh.write('\n')

        for i in range(n_rows):
            fh.write(f'{i}')
            for j in range(n_cols):
                fh.write(f'\t{mat[i][j]}')
            fh.write(f'\n')

def parse_cigar(
    s : str
) -> list[tuple[str, int]]:
    """
    parses a CIGAR string generated by A*PA aligner
    args :
        s (str) : CIGAR string
    returns :
        cigar_ops (list(tuple(str, int))) : list of (cigar_op, length) tuples 
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
    file_path : str
) -> float:
    """
    parses the alignment CSV file generated by A*PA aligner
    args :
        file_path (str) : file path to the alignment results
    returns :
        max_pident (float) : maximum percent identity
    """
    with open(file_path, 'r') as fh:
        max_pident = None
        for ln in fh:
            parts = ln.strip().split(",")
            cigar = parts[1]
            cigar_ops = parse_cigar(cigar)
            tot_l = sum([x[0] for x in cigar_ops])
            eq_l = sum([x[0] for x in cigar_ops if x[1] == '='])
            pident = eq_l / tot_l * 100
            if not max_pident:
                max_pident = pident
            else:
                if pident > max_pident:
                    max_pident = pident
    return max_pident