import os
from pathlib import Path

def _find_root(marker="pixi.toml"):
    p = Path(__file__).resolve()
    for parent in p.parents:
        if (parent / marker).exists():
            return parent
    raise RuntimeError(f"Could not find {marker} above {p}")

PROJECT_ROOT = _find_root()
DATA_ROOT = Path(os.environ.get(
    "MDPROJ_DATA_ROOT",
    "/data/bari-garnier/these/data/lbt/Bordeaux_System/Membrane_Simulation",
))
