from __future__ import annotations

from tqdm import tqdm


def iter_traj(
    universe,
    traj_slice: slice | None = None,
    desc: str | None = None,
    progress: bool = True,
):
    """Iterate over a universe's trajectory with a tqdm progress bar.

    Yields the same Timestep objects as `for ts in universe.trajectory[...]`,
    and updates the universe's atom positions as usual.
    """
    traj = universe.trajectory if traj_slice is None else universe.trajectory[traj_slice]
    return tqdm(traj, total=len(traj), desc=desc, unit="frame", disable=not progress)
