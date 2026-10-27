from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations as trans
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.lib.distances import minimize_vectors  # MDAnalysis >= 2.1
from tqdm.auto import tqdm

from .mda2cif import write_cif  # the gemmi writer from before


class KeepChainsTogether:
    """Shift whole (already unwrapped) chains by box vectors so that each one
    sits in the periodic image closest to the first chain."""

    def __init__(self, chains):
        self.chains = list(chains)

    def __call__(self, ts):
        if ts.dimensions is None or len(self.chains) < 2:
            return ts
        ref = self.chains[0].center_of_geometry()
        for ch in self.chains[1:]:
            d = ch.center_of_geometry() - ref
            d_min = minimize_vectors(d[None, :].astype(np.float32), ts.dimensions)[0]
            ch.positions = ch.positions + (d_min - d)
        return ts


def align_trajectory(
    system: Universe | AtomGroup,
    align_selection: str = "protein",
    output: str | Path | None = None,
    traj_slice: slice | None = None,
    wrap_others: bool = True,
) -> AtomGroup:
    """Make molecules whole, keep chains together, center, align. Optionally save."""
    ag_in = system.atoms if isinstance(system, Universe) else system
    u = ag_in.universe.copy()
    ag = u.atoms[ag_in.indices]
    protein = ag.select_atoms(align_selection)
    others = ag - protein
    chains = [s.atoms & protein for s in protein.segments]
    chains = [c for c in chains if len(c)]

    print(f"Aligning on '{align_selection}' ({protein.n_atoms} atoms), "
          f"{len(chains)} chain(s), {others.n_atoms} other atoms")

    # 1. TRANSFORMS ---------------------------------------------------------
    u.trajectory[0]
    pre = [trans.unwrap(ag)]                       # every bonded molecule whole
    if len(chains) > 1:
        pre.append(KeepChainsTogether(chains))     # chains in the same image
    pre.append(trans.center_in_box(protein, center="geometry"))
    if wrap_others and len(others):
        pre.append(trans.wrap(others, compound="fragments"))  # keeps molecules intact

    for t in pre:                                  # process frame 0 -> reference for fit
        t(u.trajectory.ts)
    fit = trans.fit_rot_trans(protein, protein, weights="mass")
    u.trajectory.add_transformations(*pre, fit)

    # 2. SAVE ---------------------------------------------------------------
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)

        sel = traj_slice if traj_slice is not None else slice(None)
        frame_ids = range(u.trajectory.n_frames)[sel]
        n_frames = len(frame_ids)

        # topology written from the FIRST frame that goes into the xtc
        u.trajectory[frame_ids[0]]
        cif_path = output.with_suffix(".cif")
        write_cif(ag, str(cif_path))
        print(f"Wrote topology to: {cif_path}")

        if n_frames > 1:
            xtc_path = output.with_suffix(".xtc")
            has_bonds = hasattr(ag, "bonds") and len(ag.bonds) > 0
            if has_bonds:
                b1, b2 = ag.bonds.atom1, ag.bonds.atom2
            max_bond, rgs = 0.0, []

            with mda.Writer(str(xtc_path), ag.n_atoms) as W:
                for _ in tqdm(u.trajectory[sel], total=n_frames,
                              desc=f"Writing {xtc_path.name}", unit="frame"):
                    W.write(ag)
                    if has_bonds:
                        d = np.linalg.norm(b1.positions - b2.positions, axis=1).max()
                        max_bond = max(max_bond, d)
                    rgs.append(protein.radius_of_gyration())

            print(f"Wrote aligned trajectory ({n_frames} frames) to: {xtc_path}")
            if has_bonds:
                print(f"  longest bond over all frames: {max_bond:.1f} Å (should be a few Å)")
            print(f"  protein Rg range: {min(rgs):.1f}-{max(rgs):.1f} Å "
                  f"(sudden jumps = chains drifting apart)")
        else:
            print(f"Only {n_frames} frame(s) selected, skipping .xtc output.")

    return ag
