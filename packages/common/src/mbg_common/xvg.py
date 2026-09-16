"""XVG file parsing utilities for GROMACS output."""
from pathlib import Path
import numpy as np


def read_xvg(filepath):
    """
    Parse a GROMACS .xvg file, extracting data and axis/legend metadata.

    Returns:
        dict with keys:
            'data': numpy array (n_rows, n_cols) — column 0 is typically time/step
            'labels': legend labels, one per data series
            'xlabel', 'ylabel': axis labels
            'title': plot title, or None
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"XVG file not found: {filepath}")

    data_lines = []
    labels = []
    xlabel = None
    ylabel = None
    title = None

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("@"):
                if "xaxis label" in line:
                    xlabel = line.split('"')[1] if '"' in line else None
                elif "yaxis label" in line:
                    ylabel = line.split('"')[1] if '"' in line else None
                elif "title" in line:
                    title = line.split('"')[1] if '"' in line else None
                elif line.startswith("@ s") and "legend" in line:
                    labels.append(line.split('"')[1])
                continue

            if line.startswith("#"):
                continue

            try:
                values = [float(x) for x in line.split()]
                data_lines.append(values)
            except ValueError:
                continue

    if not data_lines:
        raise ValueError(f"No data found in {filepath}")

    data = np.array(data_lines)

    if not labels:
        labels = [f"Column {i}" for i in range(data.shape[1] - 1)]

    return {
        "data": data,
        "labels": labels,
        "xlabel": xlabel or "Time",
        "ylabel": ylabel or "Value",
        "title": title,
    }
