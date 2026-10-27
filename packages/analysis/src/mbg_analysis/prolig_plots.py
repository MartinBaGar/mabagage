from pathlib import Path
import matplotlib.pyplot as plt


def plot_contact_ranges(ranges, stride=1, dt_ps=None, xlim=None, filename=None):
    rows = sorted(ranges)

    protein_names = sorted(
        {protein_name for ligand_name, protein_name in rows}
    )
    colors = plt.get_cmap("tab10")
    protein_color = {
        name: colors(i % 10)
        for i, name in enumerate(protein_names)
    }

    fig, ax = plt.subplots(figsize=(10, 0.4 * len(rows) + 1.5))

    for row_number, (ligand_name, protein_name) in enumerate(rows):
        bars = []

        for start, end, n_frames in ranges[(ligand_name, protein_name)]:
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

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{lig} - {prot}" for lig, prot in rows])
    ax.invert_yaxis()

    ax.set_xlabel("Time (ns)" if dt_ps is not None else "Frame")
    ax.set_title("Ligand-protein contacts")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=protein_color[name])
        for name in protein_names
    ]
    ax.legend(handles, protein_names, title="Protein", loc="upper right")

    if xlim is not None:
        ax.set_xlim(xlim)

    fig.tight_layout()

    if filename is not None:
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200)

    return fig
