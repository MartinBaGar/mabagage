from pathlib import Path
import matplotlib.pyplot as plt


def plot_contact_ranges(
        ranges,
        stride=1,
        dt_ps=None,
        xlim=None,
        protein_colors=None,
        filename=None,
):
    """
    Timeline plot: one row per (ligand, protein) pair, one bar per range.

    Parameters
    ----------
    ranges : dict
        Dictionary from contact_ranges().
    stride : int
        Stride used in ligands_contacts().
    dt_ps : float, optional
        Time between trajectory frames in ps. If provided, the x-axis
        is converted to ns.
    xlim : tuple, optional
        Limits in ns when dt_ps is provided, otherwise in frames.
    protein_colors : dict, optional
        Mapping from protein/chain ID to a Matplotlib color, for example:
        {"D": "tab:blue", "E": "#E69F00"}.
        Missing proteins receive colors from the tab10 colormap.
    filename : str or Path, optional
        Output image filename.
    """
    rows = sorted(ranges)
    protein_names = sorted({protein_name for _, protein_name in rows})

    # One colour per protein, so the same protein always looks the same
    protein_names = sorted({protein_name for ligand_name, protein_name in rows})
    colors = plt.get_cmap("tab10")
    
        # User-defined colors, with automatic fallback for missing chains.
    protein_colors = {} if protein_colors is None else dict(protein_colors)
    fallback = plt.get_cmap("tab10")


    protein_color = {
        name: protein_colors.get(name, fallback(i % 10))
        for i, name in enumerate(protein_names)
    }

    fig, ax = plt.subplots(
        figsize=(10, max(2.0, 0.4 * len(rows) + 1.5))
    )

    for row_number, (ligand_name, protein_name) in enumerate(rows):
        bars = []

        for start, end, n_frames in ranges[(ligand_name, protein_name)]:
            # Include the final sampled interval so one-frame contacts
            # still have visible width.
            width = end - start + stride

            if dt_ps is not None:
                start = start * dt_ps / 1000
                width = width * dt_ps / 1000

            bars.append((start, width))

        ax.broken_barh(
            bars,
            (row_number - 0.4, 0.8),
            color=protein_color[protein_name],
        )

    ax.invert_yaxis()
    
    ax.set_yticks([])
    ax.set_yticklabels("")

    # ax.set_xlabel("Time (ns)" if dt_ps is not None else "Frame")
    ax.set_xlabel("Frame", size=16)
    ax.set_ylabel(f"{ligand_name} in contact", size=16)
    ax.set_title(f"{ligand_name}-protein contacts")

    handles = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            color=protein_color[name],
            label=name,
        )
        for name in protein_names
    ]

    ax.legend(
        bbox_to_anchor=(1.1, 0.75),
        handles=handles,
        title="Ubi",
        loc="upper right",
    )

    if xlim is not None:
        ax.set_xlim(xlim)

    fig.tight_layout()

    if filename is not None:
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(filename, dpi=200, bbox_inches="tight", transparent=True)

    return fig
