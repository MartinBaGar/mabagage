import os
import re
from pymol import cmd, util
from pymol.cgo import *

# --- CONFIGURATION ---
# Add your lipid residue names and their desired color functions here
LIPID_MAP = {
    "PIPC": "cbag",  # Carbon By Atom Green
    "PWPC": "cbac",  # Cyan
    "DOPS": "cbam",  # Magenta
    "DSPS": "cbam",  # Magenta
    "IXPS": "cbam",  # Magenta
    "OEPS": "cbam",  # Magenta
    "WXPS": "cbay",  # Yellow
    "SITO": "cbas",  # Salmon
    "SXGA": "cbaw",  # White
    "SPGA": "cbaw",  # Blue
    "SSGA": "cbaw",  # Update
    "CER1": "cbaw",  # White
    "CER2": "cbab",  # Blue
}

LIPID_BONDS = {
    "GIA0": [("B0", "A1"), ("B1", "A2")],
    "GIA2": [("B0", "A1"), ("B1", "A2")],
    "GIB0": [("B0", "A1"), ("B1", "A2"), ("B2", "A3")],
}

# Define RGB for maximum precision
cmd.set_color("cb_vermillion", [0.8, 0.4, 0.0])
cmd.set_color("cb_skyblue", [0.35, 0.7, 0.9])
cmd.set_color("cb_bluegreen", [0.0, 0.6, 0.5])
cmd.set_color("cb_amber", [0.9, 0.6, 0.0])


# --- INTERNAL HELPERS (Private) ---
def _setup_spheres(selection):
    """Internal helper to standardize sphere scales for lipids."""
    cmd.show("spheres", selection)
    cmd.set("sphere_scale", 0.7, f"({selection}) and name H*")
    cmd.set("sphere_scale", 0.9, f"({selection}) and not name H*")
    
    
def _gaussian_base_name(selection, name=None):
    """
    Return an explicit name or derive one from a simple selection.

    Simple object/selection names may contain letters, numbers,
    underscores, dots, hyphens, and plus signs.
    """
    if name is not None and str(name).strip():
        base_name = str(name).strip()
    else:
        selection_name = str(selection).strip()

        if not re.fullmatch(r"[A-Za-z0-9_.+\-]+", selection_name):
            raise ValueError(
                "Cannot derive an object name from the complex selection "
                f"{selection!r}. Provide an explicit name, for example: "
                "gaussian_surface 'prot and polymer.protein', name=protein"
            )

        base_name = selection_name

    # Ensure that the generated PyMOL names are safe.
    base_name = re.sub(r"[^A-Za-z0-9_]+", "_", base_name)
    base_name = base_name.strip("_")

    if not base_name:
        raise ValueError("The generated base name is empty.")

    return base_name


def _gaussian_safe_label(value, fallback):
    """Convert a chain or segment identifier into a safe name fragment."""
    value = str(value).strip()

    if not value:
        return fallback

    value = re.sub(r"[^A-Za-z0-9_]+", "_", value)
    value = value.strip("_")

    return value or fallback


def _gaussian_identifiers(selection, split_by):
    """
    Return the unique chain or segment identifiers in a selection.
    """
    identifiers = set()

    if split_by == "chain":
        cmd.iterate(
            selection,
            "identifiers.add(chain)",
            space={"identifiers": identifiers},
        )

    elif split_by == "segi":
        cmd.iterate(
            selection,
            "identifiers.add(segi)",
            space={"identifiers": identifiers},
        )

    else:
        raise ValueError(
            "split_by must be None, 'chain', or 'segi'."
        )

    # Put blank identifiers last.
    return sorted(
        identifiers,
        key=lambda value: (value == "", value),
    )

def _copy_uniform_color(source_selection, target_object):
    """
    Copy a source selection's color to a target object if all source
    atoms share one color.

    Returns True if a color was copied, otherwise False.
    """
    color_indices = set()

    cmd.iterate(
        source_selection,
        "color_indices.add(color)",
        space={"color_indices": color_indices},
    )

    if len(color_indices) != 1:
        return False

    color_index = next(iter(color_indices))
    rgb = cmd.get_color_tuple(color_index)

    generated_color_name = f"_{target_object}_color"
    cmd.set_color(generated_color_name, rgb)
    cmd.color(generated_color_name, target_object)

    return True


# --- PUBLIC COMMANDS ---
@cmd.extend
def gaussian_surface(
    selection="all",
    name=None,
    split_by=None,
    b_floor=40.0,
    resolution=5.0,
    spacing=0.5,
    buffer=10.0,
    level=1.0,
):
    """
    Create one or more Gaussian isosurfaces from a selection.

    Parameters
    ----------
    selection : str
        Source PyMOL atom selection.

    name : str or None
        Base name for generated objects. If omitted, the name is
        derived from a simple selection or object name. Complex
        selections require an explicit name.

    split_by : None, "chain", or "segi"
        None:
            Create one Gaussian surface for the complete selection.

        "chain":
            Create one independent surface per chain identifier.

        "segi":
            Create one independent surface per segment identifier.

    b_floor : float
        Minimum B-factor used during Gaussian-map generation.
        Higher values generally create broader density.

    resolution : float
        Gaussian-map resolution in Angstroms. Higher values produce
        smoother and less detailed density.

    spacing : float
        Map grid spacing in Angstroms. Smaller values produce smoother
        surface geometry but use more memory.

    buffer : float
        Extra map-box distance around each source selection, in
        Angstroms.

    level : float
        Initial isosurface density threshold. Lower values create
        broader surfaces; higher values create tighter surfaces.

    Generated names
    ---------------
    Combined mode:
        <name>_gaussian
        <name>_gaussian_map

    Split mode:
        <name>_<identifier>_gaussian
        <name>_<identifier>_gaussian_map

    Groups:
        <name>_gaussians
        <name>_gaussian_maps

    Examples
    --------
    PyMOL command line:

        gaussian_surface prot

        gaussian_surface prot, split_by=chain

        gaussian_surface prot, name=protein, split_by=segi

        gaussian_surface "prot and polymer.protein", name=protein

        gaussian_surface vesicle, resolution=5, spacing=0.5, \
            buffer=15, level=1.0

    Python:

        gaussian_surface(
            selection="all",
            name=None,
            split_by=None,
            b_floor=40.0,
            resolution=5.0,
            spacing=1.0,
            buffer=8.0,
            level=1.0,
        )

    Notes
    -----
    The function temporarily changes gaussian_b_floor and
    gaussian_resolution, then restores their previous values.

    The generated map objects are retained so that the surface level
    can be changed afterward, for example:

        cmd.isolevel("prot_A_gaussian", 1.2)
    """

    selection = str(selection).strip()

    if not selection:
        raise ValueError("selection cannot be empty")

    if cmd.count_atoms(selection) == 0:
        raise ValueError(
            f"No atoms found in selection: {selection!r}"
        )

    base_name = _gaussian_base_name(selection, name)

    # PyMOL command-line arguments normally arrive as strings.
    b_floor = float(b_floor)
    resolution = float(resolution)
    spacing = float(spacing)
    buffer = float(buffer)
    level = float(level)

    if spacing <= 0:
        raise ValueError("spacing must be greater than zero")

    if buffer < 0:
        raise ValueError("buffer cannot be negative")

    if resolution <= 0:
        raise ValueError("resolution must be greater than zero")

    # Normalize split_by values from the PyMOL command line.
    if split_by is None:
        split_by = None
    else:
        split_by = str(split_by).strip().lower()

        if split_by in {"", "none", "false", "off", "0"}:
            split_by = None
        elif split_by not in {"chain", "segi"}:
            raise ValueError(
                "split_by must be None, 'chain', or 'segi'"
            )

    # Build the list of source selections.
    # Each item is (identifier label, PyMOL selection).
    if split_by is None:
        sources = [(None, selection)]

    else:
        sources = []

        identifiers = _gaussian_identifiers(
            selection,
            split_by,
        )

        for index, identifier in enumerate(identifiers, start=1):
            if identifier:
                source_selection = (
                    f"({selection}) and {split_by} "
                    f"{identifier}"
                )
                label = _gaussian_safe_label(
                    identifier,
                    fallback=f"{split_by}_{index}",
                )
            else:
                # Quoting an empty chain/segi can be awkward across
                # PyMOL versions. The temporary named selection makes
                # the expression unambiguous.
                source_selection = (
                    f"({selection}) and {split_by} ''"
                )
                label = f"blank_{split_by}"

            if cmd.count_atoms(source_selection) > 0:
                sources.append((label, source_selection))

        if not sources:
            raise ValueError(
                f"No valid {split_by} identifiers were found in "
                f"selection {selection!r}"
            )

    surface_group = f"{base_name}_gaussians"
    map_group = f"{base_name}_gaussian_maps"

    created_surfaces = []
    created_maps = []

    # Store the current global Gaussian settings.
    old_b_floor = cmd.get("gaussian_b_floor")
    old_resolution = cmd.get("gaussian_resolution")

    try:
        cmd.set("gaussian_b_floor", b_floor)
        cmd.set("gaussian_resolution", resolution)

        for label, source_selection in sources:
            if label is None:
                surface_name = f"{base_name}_gaussian"
                map_name = f"{base_name}_gaussian_map"
            else:
                surface_name = (
                    f"{base_name}_{label}_gaussian"
                )
                map_name = (
                    f"{base_name}_{label}_gaussian_map"
                )

            # Repeated calls replace generated objects with the same names.
            cmd.delete(surface_name)
            cmd.delete(map_name)

            cmd.map_new(
                map_name,
                "gaussian",
                spacing,
                source_selection,
                buffer,
            )

            cmd.isosurface(
                surface_name,
                map_name,
                level,
            )
            color_copied = _copy_uniform_color(
                source_selection,
                surface_name,
            )

            if not color_copied:
                print(
                    f" Note: {source_selection!r} contains multiple colors; "
                    f"{surface_name} keeps its default color."
                )

            created_surfaces.append(surface_name)
            created_maps.append(map_name)

    finally:
        # Restore the settings even if map generation raises an error.
        cmd.set("gaussian_b_floor", old_b_floor)
        cmd.set("gaussian_resolution", old_resolution)

    # Replace previous groups of the same name.
    cmd.delete(surface_group)
    cmd.delete(map_group)

    cmd.group(
        surface_group,
        " ".join(created_surfaces),
    )

    cmd.group(
        map_group,
        " ".join(created_maps),
    )

    # Retain the maps but keep them disabled.
    cmd.disable(map_group)
    cmd.rebuild()

    return created_surfaces


@cmd.extend
def color_lipids(selection="not solvent and not polymer.protein"):
    """
    Colors the bilayer by residue type using the LIPID_MAP palette.
    Usage: color_lipids [selection]
    """
    _setup_spheres(selection)

    found_any = False
    for resn, color_func in LIPID_MAP.items():
        target = f"({selection}) and resn {resn}"
        if cmd.count_atoms(target) > 0:
            # getattr lets us call util.cbag() using the string "cbag"
            getattr(util, color_func)(target, _self=cmd)
            found_any = True

    if not found_any:
        print(" Warning: No recognized lipids found in selection.")
    else:
        print(" Lipids colored successfully.")


@cmd.extend
def highlight_insertion(
    residues="133+195+197-198", protein_sel="polymer.protein"
):
    """
    Focuses on the protein-membrane interface.
    Usage: highlight_insertion [residues], [protein_selection]
    """
    # 1. Clean up view
    _setup_spheres(f"not solvent and not {protein_sel}")

    # 2. Protein representation
    cmd.dss(protein_sel)
    cmd.show_as("cartoon", protein_sel)

    # 3. Highlight specific residues
    sel_name = "insertion_site"
    cmd.select(sel_name, f"({protein_sel}) and resi {residues}")
    cmd.show("sticks", sel_name)
    cmd.color("magentas", sel_name)

    # 4. Focus
    cmd.zoom(sel_name, buffer=8)
    print(f" Focused on residues {residues} in {protein_sel}.")


@cmd.extend
def add_axes():
    w = 0.06  # cylinder width
    l = 0.75  # cylinder length
    h = 0.25  # cone hight
    d = w * 1.618  # cone base diameter

    obj = [
        CYLINDER,
        0.0,
        0.0,
        0.0,
        l,
        0.0,
        0.0,
        w,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        CYLINDER,
        0.0,
        0.0,
        0.0,
        0.0,
        l,
        0.0,
        w,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        CYLINDER,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        l,
        w,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        CONE,
        l,
        0.0,
        0.0,
        h + l,
        0.0,
        0.0,
        d,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        1.0,
        CONE,
        0.0,
        l,
        0.0,
        0.0,
        h + l,
        0.0,
        d,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        1.0,
        1.0,
        CONE,
        0.0,
        0.0,
        l,
        0.0,
        0.0,
        h + l,
        d,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        1.0,
        1.0,
    ]

    cmd.load_cgo(obj, "axes")


@cmd.extend
def vid_med():
    cmd.hide("everything", "resn W resn ION")
    pic_med()
    # --- Sphere Colors (High Contrast) ---
    cmd.set("sphere_color", "cb_vermillion", "e. A*")
    cmd.set("sphere_color", "cb_skyblue", "e. B*")
    cmd.set("sphere_color", "cb_bluegreen", "e. C*")
    cmd.set("sphere_color", "cb_amber", "e. D*")
    cmd.set("sphere_color", "yelloworange", "e. PO4")
    cmd.set("sphere_color", "black", "e. VS*")


@cmd.extend
def prep_for_movie():
    pic_med()
    # --- Sphere Colors (High Contrast) ---
    cmd.hide("sticks")
    cmd.set("sphere_color", "cb_vermillion", "e. A*")
    cmd.set("sphere_color", "cb_skyblue", "e. B*")
    cmd.set("sphere_color", "cb_bluegreen", "e. C*")
    cmd.set("sphere_color", "cb_amber", "e. D*")
    cmd.set("sphere_color", "yelloworange", "e. PO4")
    cmd.set("sphere_color", "black", "e. VS*")

    cmd.set("sphere_transparency", 0.8, "resn W resn ION")
    cmd.set("sphere_scale", 0.3, "resn W resn ION")


@cmd.extend
def pic_pretty():
    """Quick matte render for presentations."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("antialias", 1)

    # Matte settings
    cmd.set("ambient", 0.4)
    cmd.set("direct", 0.2)
    cmd.set("specular", 0)
    cmd.set("ray_shadow", 0)
    cmd.set("ray_trace_mode", 0)
    cmd.ray()


@cmd.extend
def pic_med():
    """Matte finish with depth/shadows, no outlines."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("depth_cue", 0)
    cmd.set("antialias", 2)

    # Matte settings
    cmd.set("ambient", 0.5)
    cmd.set("direct", 0.5)
    cmd.set("specular", 0.1)
    cmd.set("spec_reflect", 0)
    cmd.set("shininess", 0)
    
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadow", 1)
    cmd.set("ray_trace_mode", 0)


@cmd.extend
def pic_pub():
    """Professional Matte: High-quality geometry with thin outlines."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("antialias", 2)

    # Ultra-Matte Surface Math
    cmd.set("ambient", 0.5)
    cmd.set("direct", 0.45)
    cmd.set("specular", 0.05)
    cmd.set("spec_reflect", 0)
    cmd.set("shininess", 0)

    # High-End Geometry
    cmd.set("cartoon_sampling", 14)
    cmd.set("ribbon_sampling", 14)

    # Outlines that won't drown out small details
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.2)
    cmd.set("ray_trace_depth_factor", 0.1)
    cmd.set("ray_shadow", 1)

    cmd.ray(2400, 2400)


@cmd.extend
def show_beads():
    """Display labels indicating the bead names"""
    cmd.label("e. A*+B*+C*+VS*+D*", "'%s' % (name)")
    cmd.set("label_position", (0, 0, 20))
    cmd.set("label_size", 25)
    cmd.set("label_color", "black")
    cmd.set("sphere_transparency", 0.6)
    cmd.set("sphere_scale", 0.4)


@cmd.extend
def show_bonds():
    """Display labels indicating the bead names"""
    for lipid, pairs in LIPID_BONDS.items():
        for pair in pairs:
            bond_name = f"{lipid}_{pair[0]}-{pair[1]}"
            cmd.distance(
                bond_name,
                f"{lipid} and name {pair[0]}",
                f"{lipid} and name {pair[1]}",
            )
        cmd.group(f"{lipid}_group", f"{lipid}*")

    cmd.set("dash_gap", 0)
    cmd.set("dash_round_ends", 0)
    cmd.set("dash_width", 10)
    cmd.hide("labels")
    cmd.label("e. A*+B*+C*+VS*+D*", "'%s' % (name)")

    # --- Sphere Colors (High Contrast) ---
    cmd.set("sphere_color", "cb_vermillion", "e. A*")
    cmd.set("sphere_color", "cb_skyblue", "e. B*")
    cmd.set("sphere_color", "cb_bluegreen", "e. C*")
    cmd.set("sphere_color", "cb_amber", "e. D*")
    cmd.set("sphere_color", "black", "e. VS*")

    # --- Label Colors (Maximum Contrast vs Sphere) ---
    # Dark text for light spheres, White text for dark spheres
    cmd.set("label_color", "black", "e. A*")  # Black on Vermillion
    cmd.set("label_color", "black", "e. B*")  # Black on Sky Blue
    cmd.set("label_color", "white", "e. C*")  # White on Blue-Green
    cmd.set("label_color", "black", "e. D*")  # Black on Amber
    cmd.set("label_color", "white", "e. VS*")  # White on Black

    # --- Global Settings ---
    cmd.set("label_position", (0, 0, 10))
    cmd.set("label_size", 20)
    cmd.set("label_font_id", 7)  # Bold font for better readability
    cmd.set("sphere_transparency", 0.6)
    cmd.set("sphere_scale", 0.3)
    cmd.set("sphere_scale", 0.2, "e. VS*")


@cmd.extend
def my_help():
    """
    Displays this help menu.
    """
    print("\n" + "=" * 50)
    print(" LAB VISUALIZATION LIBRARY ".center(50, " "))
    print("=" * 50)

    # List of functions to document
    commands = [color_lipids, highlight_insertion, my_help, reload_lib]

    for func in commands:
        name = func.__name__.ljust(20)
        # Use the first line of the docstring as the description
        desc = (
            func.__doc__.strip().split("\n")[0]
            if func.__doc__
            else "No description."
        )
        print(f" {name} : {desc}")

    print("=" * 50 + "\n")


@cmd.extend
def reload_lib():
    """
    Reloads this script after you make changes to the code.
    """
    script_path = os.path.abspath(__file__)
    cmd.run(script_path)
    print(f" Reloaded: {script_path}")


# --- INITIALIZATION ---
# This runs only when the file is first loaded or reloaded
# if __name__ == "__main__" or "pymol" in globals():
#     my_help()

#     print
