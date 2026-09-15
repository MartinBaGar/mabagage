"""Plot GROMACS .xvg files with automatic label detection."""
from pathlib import Path
import argparse
import sys

import matplotlib.pyplot as plt

from mdproj_common.xvg import read_xvg
from mdproj_common.plotting import apply_md_style


def plot_xvg(xvg_file, output=None, column=1, **kwargs):
    """
    Plot a single data column from a GROMACS .xvg file.

    Args:
        xvg_file: path to .xvg file
        output: output filename (auto-generated from input if None)
        column: which data column to plot (default 1, after time/step column 0)
        **kwargs: optional overrides — xlabel, ylabel, title, xlim, ylim, figsize

    Returns:
        Path to the saved figure
    """
    data_dict = read_xvg(xvg_file)
    data = data_dict["data"]

    if output is None:
        output = Path(xvg_file).stem + "_plot.png"

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    apply_md_style()
    figsize = kwargs.get("figsize", (8, 5))
    fig, ax = plt.subplots(figsize=figsize)

    time = data[:, 0]
    values = data[:, column]
    ax.plot(time, values, linewidth=1.5)

    ax.set_xlabel(kwargs.get("xlabel") or data_dict["xlabel"])
    ax.set_ylabel(kwargs.get("ylabel") or data_dict["ylabel"])

    if kwargs.get("title"):
        ax.set_title(kwargs["title"])
    if kwargs.get("xlim"):
        ax.set_xlim(kwargs["xlim"])
    if kwargs.get("ylim"):
        ax.set_ylim(kwargs["ylim"])

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot GROMACS .xvg files with automatic label detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m mdproj_analysis.xvg_plots energy.xvg
  python -m mdproj_analysis.xvg_plots energy.xvg -o my_energy.png
  python -m mdproj_analysis.xvg_plots pressure.xvg --ylabel "Pressure (bar)"
  python -m mdproj_analysis.xvg_plots temperature.xvg --xlim 0 1000
  python -m mdproj_analysis.xvg_plots energy.xvg --title "NVT Equilibration"
        """,
    )
    parser.add_argument("xvg_file", help="Input .xvg file")
    parser.add_argument("-o", "--output", help="Output filename (auto-generated if not specified)")
    parser.add_argument("--column", type=int, default=1, help="Data column to plot (default: 1)")
    parser.add_argument("--xlabel", help="Override X-axis label")
    parser.add_argument("--ylabel", help="Override Y-axis label")
    parser.add_argument("--title", help="Add plot title")
    parser.add_argument("--xlim", type=float, nargs=2, metavar=("MIN", "MAX"), help="X-axis limits")
    parser.add_argument("--ylim", type=float, nargs=2, metavar=("MIN", "MAX"), help="Y-axis limits")
    parser.add_argument("--figsize", type=float, nargs=2, default=(8, 5), metavar=("WIDTH", "HEIGHT"), help="Figure size in inches")

    args = parser.parse_args()

    plot_kwargs = {
        "column": args.column,
        "xlabel": args.xlabel,
        "ylabel": args.ylabel,
        "title": args.title,
        "xlim": args.xlim,
        "ylim": args.ylim,
        "figsize": tuple(args.figsize),
    }
    plot_kwargs = {k: v for k, v in plot_kwargs.items() if v is not None}

    try:
        result = plot_xvg(args.xvg_file, args.output, **plot_kwargs)
        print(result, end="")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
