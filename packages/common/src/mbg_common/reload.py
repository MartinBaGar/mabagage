import sys
import importlib
import os

_mtimes = {}


def reload_mbg(prefix="mbg_", verbose=True):
    """Reload all currently-imported modules under `prefix`. Always
    reloads matched modules (cheap); only prints when source actually
    changed since the last call, so it's safe to call before every
    block (e.g. via org's :prologue).
    """
    mods = [
        name for name in sys.modules if name.split(".")[0].startswith(prefix)
    ]
    mods.sort(key=lambda n: n.count("."))

    changed = []
    for name in mods:
        mod = sys.modules.get(name)
        path = getattr(mod, "__file__", None)
        if not path:
            continue
        try:
            mtime = os.stat(path).st_mtime
        except OSError:
            continue
        prev = _mtimes.get(path)
        if prev is not None and prev != mtime:
            changed.append(name)
        _mtimes[path] = mtime

    for _ in range(2):
        for name in mods:  # reload everything, not just `changed`
            mod = sys.modules.get(name)
            if mod is None:
                continue
            try:
                importlib.reload(mod)
            except Exception as e:
                if verbose:
                    print(f"[reload_mbg] failed on {name}: {e}")

    if verbose and changed:
        print(
            f"[reload_mbg] reloaded {len(changed)} module(s): {', '.join(changed)}"
        )
