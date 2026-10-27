"""Export a solvent-free Martini snapshot to mmCIF using Gemmi."""

import warnings

import gemmi
import numpy as np


def write_cif(atomgroup, path):
    """Write the current frame; return global atom indices in file order.

    Assumes:
      - chainIDs are already assigned.
      - Polymer residue names are recognized by MDAnalysis.
      - Selected polymer residues define the complete sequence.
      - Each lipid/small molecule occupies one residue.
      - Residue identifiers are unique within each author chain
        (checked: a ValueError is raised otherwise).
      - No solvent is included.

    Preserves chainIDs, resids, residue names, and bead names.
    Does not export Martini masses, force-field types, or bonds.

    Raises ValueError for empty selections, missing/empty chain IDs,
    empty residue or bead names, non-finite coordinates, and duplicate
    label_asym_id values.
    """
    ag = atomgroup.atoms

    # ---- Input validation ------------------------------------------------
    if ag.n_atoms == 0:
        raise ValueError("Cannot export an empty AtomGroup.")

    try:
        chain_ids = ag.chainIDs
    except AttributeError as exc:
        raise ValueError(
            "The topology has no 'chainIDs' attribute; assign chain IDs "
            "before exporting."
        ) from exc

    names = ag.names
    resnames = ag.resnames
    resids = ag.resids
    positions = ag.positions

    if any(not str(c).strip() for c in np.unique(chain_ids.astype(str))):
        n_bad = int(np.sum(np.char.strip(chain_ids.astype(str)) == ""))
        raise ValueError(f"{n_bad} selected atoms have an empty chain ID.")

    if any(not str(r).strip() for r in np.unique(resnames.astype(str))):
        raise ValueError("Some selected atoms have an empty residue name.")

    if any(not str(n).strip() for n in np.unique(names.astype(str))):
        raise ValueError("Some selected atoms have an empty bead name.")

    if not np.all(np.isfinite(positions)):
        raise ValueError("Some coordinates are NaN or infinite.")

    # ---- Polymer classification -------------------------------------------
    is_polymer = np.isin(
        ag.indices,
        ag.select_atoms("protein or nucleic").indices,
    )

    # ---- Group selected beads by their original residue -------------------
    groups = {}
    for i, resindex in enumerate(ag.resindices):
        groups.setdefault(int(resindex), []).append(i)

    chains = {}
    written = []

    for indices in groups.values():
        first = indices[0]
        chain_id = str(chain_ids[first])

        residue = gemmi.Residue()
        residue.name = str(resnames[first])
        residue.seqid = gemmi.SeqId(int(resids[first]), " ")
        residue.entity_type = (
            gemmi.EntityType.Polymer
            if is_polymer[first]
            else gemmi.EntityType.NonPolymer
        )

        for i in indices:
            bead = gemmi.Atom()
            bead.name = str(names[i])
            bead.element = gemmi.Element("X")
            bead.pos = gemmi.Position(*map(float, positions[i]))
            bead.occ = 1.0
            residue.add_atom(bead)

        chains.setdefault(chain_id, []).append((residue, indices))

    try:
        model = gemmi.Model(1)
    except TypeError:
        model = gemmi.Model("1")

    for chain_id, residues in chains.items():
        chain = gemmi.Chain(chain_id)
        for residue, indices in residues:
            chain.add_residue(residue)
            written.extend(indices)
        model.add_chain(chain)

    structure = gemmi.Structure()
    structure.name = "cg_system"
    structure.add_model(model)

    # ---- Periodic box (only if it is physically valid) --------------------
    box = ag.universe.dimensions
    if box is not None:
        box = np.asarray(box, dtype=float)
        valid_box = (
            box.shape == (6,)
            and np.all(np.isfinite(box))
            and np.all(box[:3] > 0)
            and np.all((box[3:] > 0) & (box[3:] < 180))
        )
        if valid_box:
            structure.cell = gemmi.UnitCell(*map(float, box))
            structure.spacegroup_hm = "P 1"
        else:
            warnings.warn(
                f"Invalid periodic box {box.tolist()}; "
                "writing the file without a unit cell.",
                stacklevel=2,
            )

    # ---- Entities ---------------------------------------------------------
    # Keep Gemmi's individual non-polymer subchains.
    structure.setup_entities()

    # Use the selected polymer residues as the complete sequences.
    for entity in structure.entities:
        if entity.entity_type == gemmi.EntityType.Polymer:
            span = structure[0].get_subchain(entity.subchains[0])
            entity.full_sequence = [residue.name for residue in span]

    structure.deduplicate_entities()

    # ---- Rename label instances -------------------------------------------
    # renamed: old Gemmi subchain -> new label_asym_id
    # owner:   new label_asym_id  -> old Gemmi subchain (collision check)
    renamed = {}
    owner = {}

    for chain in structure[0]:
        for residue in chain:
            old = residue.subchain

            if residue.entity_type == gemmi.EntityType.Polymer:
                new = chain.name
            else:
                new = f"{chain.name}_{residue.seqid.num}"

            if owner.setdefault(new, old) != old:
                raise ValueError(
                    f"Duplicate label_asym_id {new!r}: residue identifiers "
                    f"are not unique in chain {chain.name!r} "
                    f"(resid {residue.seqid.num}, residue {residue.name!r})."
                )

            renamed[old] = new
            residue.subchain = new

    # Keep entity references consistent with the renamed instances.
    for entity in structure.entities:
        entity.subchains = [renamed[s] for s in entity.subchains]

    # ---- Numeric entity IDs, descriptions, molecule counts ----------------
    entity_by_subchain = {}
    descriptions = {}
    molecule_counts = {}

    for number, entity in enumerate(structure.entities, start=1):
        entity.name = str(number)
        molecule_counts[entity.name] = len(entity.subchains)

        for subchain in entity.subchains:
            entity_by_subchain[subchain] = entity.name

    for chain in structure[0]:
        for residue in chain:
            entity_id = entity_by_subchain[residue.subchain]
            residue.entity_id = entity_id

            if residue.entity_type == gemmi.EntityType.Polymer:
                description = f"Protein {chain.name}"
            else:
                description = residue.name

            descriptions.setdefault(entity_id, description)

    # ---- Finalize and write -----------------------------------------------
    structure.assign_label_seq_id()
    structure.assign_het_flags()
    structure.assign_serial_numbers()

    document = structure.make_mmcif_document()
    block = document.sole_block()

    # Add metadata without replacing other generated entity fields.
    entity_table = block.get_mmcif_category("_entity.")
    ids = [str(x) for x in entity_table["id"]]
    entity_table["pdbx_description"] = [descriptions[i] for i in ids]
    entity_table["pdbx_number_of_molecules"] = [
        str(molecule_counts[i]) for i in ids
    ]
    block.set_mmcif_category("_entity.", entity_table)

    document.write_file(str(path))
    return ag.indices[np.asarray(written, dtype=np.int64)].copy()
