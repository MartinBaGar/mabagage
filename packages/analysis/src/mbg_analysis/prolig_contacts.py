from MDAnalysis.lib.distances import capped_distance

from mbg_common.iter_traj import iter_traj


def get_protein_name(residue, protein_id="segid"):
    """
    Return the name of the protein a residue belongs to.

    protein_id="segid"   -> use the segment name   (typical for GROMACS TPR files)
    protein_id="chainID" -> use the chain ID       (typical for PDB files)
    """
    if protein_id == "segid":
        return residue.segment.segid
    if protein_id == "chainID":
        return residue.atoms[0].chainID
    raise ValueError("protein_id must be 'segid' or 'chainID'")


def ligands_contacts(universe, protein, ligands, cutoff=7, stride=1, protein_id="segid"):
    """
    For every frame and every ligand, find which protein it touches
    and through which residues.

    Returns a nested dictionary:
        contacts[frame][ligand_name][protein_name] = [(resname, resid), ...]

    A ligand with no contact in a frame gets an empty dictionary.
    """
    contacts = {}

    for ts in iter_traj(universe, slice(None, None, stride), desc="Ligand contacts"):

        # 1. Start every ligand with an empty dictionary (so ligands without contact still appear)
        frame_contacts = {}
        for ligand in ligands.residues:
            frame_contacts[f"{ligand.resname}_{ligand.resindex}"] = {}

        # 2. Find all (protein atom, ligand atom) pairs closer than the cutoff
        pairs = capped_distance(
            protein.positions,
            ligands.positions,
            max_cutoff=cutoff,
            box=universe.dimensions,
            return_distances=False,
        )
        
        # 3. Turn each atom pair into "this ligand touches this residue of this protein"
        for protein_index, ligand_index in pairs:
            residue = protein[protein_index].residue
            ligand = ligands[ligand_index].residue

            ligand_name = f"{ligand.resname}_{ligand.resindex}"
            protein_name = get_protein_name(residue, protein_id)

            # setdefault creates an empty set the first time we see this protein
            residues_touched = frame_contacts[ligand_name].setdefault(protein_name, set())
            residues_touched.add((residue.resname, residue.resid))

        # 4. Store the result for this frame (sets become sorted lists)
        contacts[ts.frame] = {}
        for ligand_name, proteins in frame_contacts.items():
            contacts[ts.frame][ligand_name] = {
                protein_name: sorted(residues)
                for protein_name, residues in proteins.items()
            }

    return contacts

def contact_ranges(contacts):
    """
    Turn per-frame contacts into ranges of consecutive frames.

    Input:  the dictionary from ligands_contacts()
    Output: ranges[(ligand_name, protein_name)] = [(start_frame, end_frame, n_frames), ...]

    A range ends as soon as the ligand stops touching that protein
    in one of the analysed frames.
    """
    ranges = {}       # finished ranges
    open_ranges = {}  # ranges still running: key -> [start_frame, last_frame, n_frames]

    for frame, ligand_dict in contacts.items():

        # 1. Which (ligand, protein) pairs are touching in this frame?
        touching_now = set()
        for ligand_name, proteins_touched in ligand_dict.items():
            for protein_name in proteins_touched:
                touching_now.add((ligand_name, protein_name))

        # 2. Extend the ranges that continue, start new ones for new contacts
        for key in touching_now:
            if key in open_ranges:
                open_ranges[key][1] = frame   # new last frame
                open_ranges[key][2] += 1      # one more frame in the range
            else:
                open_ranges[key] = [frame, frame, 1]

        # 3. Close the ranges whose contact just stopped
        for key in list(open_ranges):
            if key not in touching_now:
                finished = open_ranges.pop(key)
                ranges.setdefault(key, []).append(tuple(finished))

    # 4. Close whatever is still running at the end of the trajectory
    for key, unfinished in open_ranges.items():
        ranges.setdefault(key, []).append(tuple(unfinished))

    return ranges
