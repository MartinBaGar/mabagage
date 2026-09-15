from pathlib import Path
from typing import Union, Optional, Tuple

import MDAnalysis as mda
from MDAnalysis import transformations as trans
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
from tqdm import tqdm


def align_trajectory(
    atomgroup: Union[Universe, AtomGroup],
    out_prefix: Optional[Union[str, Path]] = None,
    stride: int = 100,
) -> Optional[Tuple[Path, Path]]:
    """
    Fit-and-rotate-align a trajectory in memory.
    Optionally write the aligned trajectory and structure to disk.

    Parameters
    ----------
    atomgroup : MDAnalysis.Universe or MDAnalysis.AtomGroup
        The system to align. If an AtomGroup is provided, the alignment
        is calculated based on these atoms, but applied to the whole Universe.
    out_prefix : str or Path, optional
        The file path prefix for saving. If provided, writes out_prefix.xtc
        and out_prefix.gro. If None, only aligns in memory.
    stride : int, default=100
        Write every Nth frame if saving to disk.

    Returns
    -------
    Tuple[Path, Path] or None
        Paths to the saved (xtc, gro) files if out_prefix is given, else None.
    """
    # 1. Standardize input (handle both Universe and AtomGroup)
    if isinstance(atomgroup, Universe):
        ag = atomgroup.atoms
        u = atomgroup
    else:
        ag = atomgroup
        u = atomgroup.universe

    # 2. Reset to the first frame to use it as the reference structure
    u.trajectory[0]

    # 3. Create and apply the alignment transformation
    align_transform = trans.fit_rot_trans(ag, ag, weights=ag.masses)
    u.trajectory.add_transformations(align_transform)

    # 4. If no save path is provided, we are done
    if out_prefix is None:
        print("Trajectory aligned in memory (no files saved).")
        return None

    # 5. Handle writing to disk
    out_prefix = Path(out_prefix)

    # Ensure the target directory exists
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    output_xtc = out_prefix.with_suffix(".xtc")
    output_gro = out_prefix.with_suffix(".gro")

    print(f"Writing trajectory to {output_xtc}...")
    with mda.Writer(str(output_xtc), ag.n_atoms) as W:
        for _ in tqdm(
            u.trajectory[::stride], desc="Writing frames", unit="frames"
        ):
            W.write(ag)

    print(f"Writing structure to {output_gro}...")
    ag.write(str(output_gro))
    print("Done!")

    return output_xtc, output_gro
