from pathlib import Path
from typing import Union, Optional

import MDAnalysis as mda
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
from tqdm import tqdm


def get_nosolv(
    topology: Optional[Union[str, Path]] = None,
    trajectory: Optional[str] = None,
    universe: Optional[Universe] = None,
    output: Optional[Union[str, Path]] = None,
) -> AtomGroup:
    """
    Get an AtomGroup containing no solvent or ions.

    Parameters
    ----------
    topology, trajectory, universe :
        As before — provide either a topology (+ optional trajectory) or
        an existing Universe.
    output :
        Optional path/prefix to save the solvent-free system to disk.
        If provided, writes:
          - "{output}.gro"  (solvent-free topology/first frame)
          - "{output}.xtc"  (solvent-free trajectory, only written if
            the universe has multiple frames)
        The extensions are added automatically.

    Returns
    -------
    AtomGroup
        Solvent-free AtomGroup (still linked to the original universe).
    """
    if universe is None:
        if trajectory is not None:
            universe = mda.Universe(str(topology), str(trajectory))
        else:
            universe = mda.Universe(str(topology))

    solvent_resnames = {"W", "NA", "CL", "K", "ION"}
    selection_string = (
        "not ( " + " or ".join(f"resname {r}" for r in solvent_resnames) + " )"
    )
    print(
        f"Removed solvent from the system using the following selection string: \n{selection_string}"
    )

    nosolv_ag = universe.select_atoms(selection_string)

    if output is not None:
        output = Path(output)
        gro_path = output.with_suffix(".gro")
        nosolv_ag.write(str(gro_path))
        print(f"Wrote solvent-free topology to: {gro_path}")

        n_frames = universe.trajectory.n_frames
        if n_frames > 1:
            xtc_path = output.with_suffix(".xtc")
            with mda.Writer(str(xtc_path), nosolv_ag.n_atoms) as W:
                for ts in tqdm(
                    universe.trajectory,
                    total=n_frames,
                    desc=f"Writing {xtc_path.name}",
                    unit="frame",
                ):
                    W.write(nosolv_ag)
            print(f"Wrote solvent-free trajectory ({n_frames} frames) to: {xtc_path}")
        else:
            print("Universe has only 1 frame — skipping trajectory (.xtc) output.")

    return nosolv_ag
