import matplotlib.pyplot as plt


def plot_contact_ranges(ranges, stride=1, dt_ps=None, filename=None):
    """
    Timeline plot: one row per (ligand, protein) pair, one bar per range.

    ranges   : dictionary from contact_ranges()
    stride   : the stride you used in ligand_contacts()
    dt_ps    : time between two trajectory frames in ps (u.trajectory.dt).
               If given, the x axis is in ns, otherwise in frames.
    filename : if given, the figure is saved there (e.g. "ranges.png")
    """
    rows = sorted(ranges)   # one row per (ligand_name, protein_name)

    # One colour per protein, so the same protein always looks the same
    protein_names = sorted({protein_name for ligand_name, protein_name in rows})
    colors = plt.get_cmap("tab10")
    protein_color = {name: colors(i % 10) for i, name in enumerate(protein_names)}

    fig, ax = plt.subplots(figsize=(10, 0.4 * len(rows) + 1.5))

    for row_number, (ligand_name, protein_name) in enumerate(rows):
        bars = []
        for start, end, n_frames in ranges[(ligand_name, protein_name)]:
            # +stride so that a one-frame range still has a visible width
            width = end - start + stride
            if dt_ps is not None:
                start, width = start * dt_ps / 1000, width * dt_ps / 1000  # ps -> ns
            bars.append((start, width))

        ax.broken_barh(bars, (row_number - 0.4, 0.8), color=protein_color[protein_name])

    # Row labels
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{lig} - {prot}" for lig, prot in rows])
    ax.invert_yaxis()   # first row at the top

    ax.set_xlabel("Time (ns)" if dt_ps is not None else "Frame")
    ax.set_title("Ligand-protein contacts")

    # Legend: protein name -> colour
    handles = [plt.Rectangle((0, 0), 1, 1, color=protein_color[n]) for n in protein_names]
    ax.legend(handles, protein_names, title="Protein", loc="upper right")

    fig.tight_layout()
    if filename is not None:
        fig.savefig(filename, dpi=200)
    return fig
