from pathlib import Path
from typing import Union, Optional

import MDAnalysis as mda
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup


def get_nosolv(
    topology: Optional[Union[str, Path]] = None,
    trajectory: Optional[str] = None,
    universe: Optional[Universe] = None,
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

    return universe.select_atoms(selection_string)
