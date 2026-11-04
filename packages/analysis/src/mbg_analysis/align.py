import string
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis import transformations as trans
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe
from MDAnalysis.lib.distances import minimize_vectors
from tqdm.auto import tqdm


class KeepChainsTogether:
    """Keep protein chains in the same periodic image.

    Each chain is assumed to have already been unwrapped. Whole chains are
    shifted by box vectors so that they occupy the periodic image closest to
    the first chain.
    """

    def __init__(self, chains):
        self.chains = list(chains)

    def __call__(self, ts):
        if ts.dimensions is None or len(self.chains) < 2:
            return ts

        reference = self.chains[0].center_of_geometry()

        for chain in self.chains[1:]:
            displacement = chain.center_of_geometry() - reference
            minimum_displacement = minimize_vectors(
                displacement[None, :].astype(np.float32),
                ts.dimensions,
            )[0]

            chain.positions += minimum_displacement - displacement

        return ts


def align_trajectory(
    system: Universe | AtomGroup,
    align_selection: str = "protein",
    output: str | Path | None = None,
    traj_slice: slice | None = None,
    wrap_others: bool = True,
) -> AtomGroup:
    """Unwrap, center, and align a trajectory in memory.

    Molecules are made whole, protein chains are kept together, and the
    selected atoms are centered and aligned against the first frame.

    If ``output`` is provided, a PDB topology is written. An XTC trajectory is
    also written when more than one frame is selected.

    The input Universe is left untouched because transformations are applied
    to a copy.

    Parameters
    ----------
    system
        Universe or AtomGroup to process.
    align_selection
        MDAnalysis selection used for centering and alignment.
    output
        Base output path. Its suffix is replaced by ``.pdb`` and ``.xtc``.
    traj_slice
        Optional slice selecting trajectory frames to write.
    wrap_others
        If true, non-alignment atoms are wrapped by molecular fragment.

    Returns
    -------
    AtomGroup
        AtomGroup belonging to the copied and transformed Universe.
    """
    # Work on a copy because MDAnalysis transformations can only be added
    # once to a Universe.
    ag_in = system.atoms if isinstance(system, Universe) else system
    u = ag_in.universe.copy()
    ag = u.atoms[ag_in.indices]

    align_group = ag.select_atoms(align_selection)
    if align_group.n_atoms == 0:
        raise ValueError(
            f"The alignment selection {align_selection!r} contains no atoms."
        )

    others = ag - align_group

    chains = [
        segment.atoms & align_group
        for segment in align_group.segments
    ]
    chains = [chain for chain in chains if chain.n_atoms > 0]

    print(
        f"Aligning on {align_selection!r} "
        f"({align_group.n_atoms} atoms), "
        f"{len(chains)} chain(s), "
        f"{others.n_atoms} other atoms"
    )

    # Initialize frame 0, which becomes the alignment reference.
    u.trajectory[0]

    transforms = [trans.unwrap(ag)]

    if len(chains) > 1:
        transforms.append(KeepChainsTogether(chains))

    transforms.append(
        trans.center_in_box(align_group, center="geometry")
    )

    if wrap_others and others.n_atoms > 0:
        transforms.append(
            trans.wrap(others, compound="fragments")
        )

    # Process frame 0 before constructing the fitting transformation so that
    # its coordinates represent the desired unwrapped/centered reference.
    for transformation in transforms:
        transformation(u.trajectory.ts)

    fit = trans.fit_rot_trans(
        align_group,
        align_group,
        weights="mass",
    )

    u.trajectory.add_transformations(*transforms, fit)

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)

        selected_indices = list(
            range(u.trajectory.n_frames)[
                traj_slice if traj_slice is not None else slice(None)
            ]
        )
        n_frames = len(selected_indices)

        if n_frames == 0:
            raise ValueError("The requested trajectory slice contains no frames.")

        # Write the topology using the first selected trajectory frame.
        u.trajectory[selected_indices[0]]

        # Assign chain IDs to protein segments for better PDB/PyMOL handling.
        if not hasattr(u.atoms, "chainIDs"):
            u.add_TopologyAttr("chainIDs")

        alphabet = (
            string.ascii_uppercase
            + string.ascii_lowercase
            + string.digits
        )
        protein_ag = ag.select_atoms("protein")

        for index, segment in enumerate(protein_ag.segments):
            chain_id = alphabet[index % len(alphabet)]
            target_atoms = segment.atoms & protein_ag
            target_atoms.chainIDs = chain_id
            target_atoms.segments.segids = f"{chain_id}"

        pdb_path = output.with_suffix(".pdb")
        ag.write(str(pdb_path), bonds=None)
        print(f"Wrote aligned topology to: {pdb_path}")

        if n_frames > 1:
            xtc_path = output.with_suffix(".xtc")

            with mda.Writer(str(xtc_path), ag.n_atoms) as writer:
                for frame_index in tqdm(
                    selected_indices,
                    total=n_frames,
                    desc=f"Writing {xtc_path.name}",
                    unit="frame",
                ):
                    u.trajectory[frame_index]
                    writer.write(ag)

            print(
                f"Wrote aligned trajectory ({n_frames} frames) to: "
                f"{xtc_path}"
            )
        else:
            print(
                f"Only {n_frames} frame(s) selected — "
                "skipping .xtc output."
            )

    return ag
