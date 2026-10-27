from pathlib import Path

import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
from biotite.structure.io.pdb import PDBFile

def convert_struc_to_fasta(
        struc_file: Path | str,
):
    """
    Write a fasta file from a structure
    """
    struc = PDBFile.read(struc_file)
    return struc
