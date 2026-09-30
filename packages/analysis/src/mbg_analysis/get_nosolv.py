import string
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
    traj_slice: Optional[slice] = None,
) -> AtomGroup:
    """
    Get an AtomGroup containing no solvent or ions.
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
        
        # --- CHAIN ID ASSIGNMENT  ---
        if not hasattr(universe.atoms, "chainIDs"):
            universe.add_TopologyAttr("chainIDs")
            
        alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
        
        # Use GROMACS's native segments instead of guessing bonds via fragments
        protein_ag = nosolv_ag.select_atoms("protein")
        for i, seg in enumerate(protein_ag.segments):
            # Assign the chain ID only to the protein atoms within this segment
            (seg.atoms & protein_ag).chainIDs = alphabet[i % len(alphabet)]
            
        print(f"Assigned Chain IDs across {len(protein_ag.segments)} protein segment(s).")
        # ---------------------------------------------------

        pdb_path = output.with_suffix(".pdb")
        nosolv_ag.write(str(pdb_path))
        print(f"Wrote solvent-free topology to: {pdb_path}")

        # --- SLICING LOGIC ---
        if traj_slice is not None:
            frames_to_write = universe.trajectory[traj_slice]
            n_frames = len(frames_to_write)
        else:
            frames_to_write = universe.trajectory
            n_frames = universe.trajectory.n_frames
        # ---------------------

        if n_frames > 1:
            xtc_path = output.with_suffix(".xtc")
            with mda.Writer(str(xtc_path), nosolv_ag.n_atoms) as W:
                for ts in tqdm(
                    frames_to_write,
                    total=n_frames,
                    desc=f"Writing {xtc_path.name}",
                    unit="frame",
                ):
                    W.write(nosolv_ag)
            print(f"Wrote solvent-free trajectory ({n_frames} frames) to: {xtc_path}")
        else:
            print(f"Only {n_frames} frame(s) selected — skipping trajectory (.xtc) output.")

    return nosolv_ag
