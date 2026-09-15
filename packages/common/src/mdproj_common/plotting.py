"""Shared plotting style and Org-Babel output helpers."""
import os
import matplotlib.pyplot as plt


def apply_md_style():
    """Modern style for all MD plots."""
    plt.rcParams.update(
        {
            "axes.prop_cycle": plt.cycler(color=["#1B9E77", "#D95F02", "#7570B3"]),
            "axes.grid": True,
            "grid.color": "#E5E5E5",
            "grid.linestyle": "--",
            "grid.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "figure.facecolor": "white",
        }
    )


def save_and_link(fig, path):
    """Save a figure and print an Org-mode file link — for functions that
    return an unsaved fig rather than saving internally themselves."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[[file:{path}]]")
