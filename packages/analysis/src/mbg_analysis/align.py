import string
from pathlib import Path
from typing import Union, Optional

import MDAnalysis as mda
from MDAnalysis import transformations as trans
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
from tqdm.auto import tqdm


def align_trajectory(
    system: Union[Universe, AtomGroup],
    align_selection: str = "protein",
    output: Optional[Union[str, Path]] = None,
    traj_slice: Optional[slice] = None,
) -> AtomGroup:
    """
    Unwrap PBCs and align a trajectory in-memory. Optionally save to disk.
    """
    # 1. SETUP: Standardize input to an AtomGroup and get the parent Universe
    ag = system.atoms if isinstance(system, Universe) else system
    u = ag.universe
    align_group = ag.select_atoms(align_selection)
    
    print(f"Aligning in memory based on: '{align_selection}' ({align_group.n_atoms} atoms)")

    # 2. TRANSFORM: Set up on-the-fly unwrapping and alignment
    u.trajectory[0]  # Jump to frame 0 to use as the reference structure
    
    unwrap = trans.unwrap(align_group)
    unwrap(u.trajectory.ts)  # Initialize unwrapping reference
    
    fit = trans.fit_rot_trans(align_group, align_group, weights="mass")
    
    u.trajectory.add_transformations(unwrap, fit)

    # 3. SAVE (Optional): Mirroring the get_nosolv logic
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)

        # --- PYMOL/MARTINI FIXES ---
        if not hasattr(u.atoms, "chainIDs"):
            u.add_TopologyAttr("chainIDs")
            
        alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
        protein_ag = ag.select_atoms("protein")
        
        for i, seg in enumerate(protein_ag.segments):
            chain_letter = alphabet[i % len(alphabet)]
            target_atoms = seg.atoms & protein_ag
            
            target_atoms.chainIDs = chain_letter
            target_atoms.segids = f"P{chain_letter}"  # Forces TER boundaries for PyMOL
            
        pdb_path = output.with_suffix(".pdb")
        ag.write(str(pdb_path), bonds=None)  # Prevents chaotic PyMOL bond guessing
        print(f"Wrote aligned topology to: {pdb_path}")
        # ---------------------------

        # --- SLICING & WRITING ---
        frames_to_write = u.trajectory[traj_slice] if traj_slice else u.trajectory
        n_frames = len(frames_to_write) if traj_slice else u.trajectory.n_frames

        if n_frames > 1:
            xtc_path = output.with_suffix(".xtc")
            with mda.Writer(str(xtc_path), ag.n_atoms) as W:
                for _ in tqdm(frames_to_write, total=n_frames, desc=f"Writing {xtc_path.name}", unit="frame"):
                    W.write(ag)
            print(f"Wrote aligned trajectory ({n_frames} frames) to: {xtc_path}")
        else:
            print(f"Only {n_frames} frame(s) selected — skipping .xtc output.")

    return ag
