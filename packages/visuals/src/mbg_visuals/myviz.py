import os
import re
import tempfile
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

# --- OUTLINE HELPERS (Private) ---
_OUTLINE_COLOR_KEYS = (
    "surface_color", "cartoon_color", "stick_color", "sphere_color",
    "line_color", "mesh_color", "dot_color", "ribbon_color",
)
_OUTLINE_TRANS_KEYS = (
    "transparency", "cartoon_transparency",
    "stick_transparency", "sphere_transparency",
)
_OUTLINE_RECOLORABLE = {
    "object:surface", "object:mesh", "object:slice", "object:map",
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


def _outline_force_opaque():
    """Remove transparency at atom, object and global level."""
    objects = cmd.get_names("objects")
    for key in _OUTLINE_TRANS_KEYS:
        cmd.unset(key, "(all)")
        for obj in objects:
            cmd.unset(key, obj)
        cmd.set(key, 0)


def _outline_flat_render(path, selection, white_objects, width, height,
                         keep_transparency, view):
    """
    Render a flat mask: `selection` (and `white_objects`) in white,
    everything else in black, no shading. Pixel brightness equals
    the per-pixel visibility of the selection.
    """
    if not keep_transparency:
        _outline_force_opaque()

    # Molecules
    for key in _OUTLINE_COLOR_KEYS:
        cmd.set(key, "black", "all")
        cmd.set(key, "white", selection)
    cmd.color("black", "all")
    cmd.color("white", selection)
    for key in _OUTLINE_TRANS_KEYS:
        cmd.set(key, 0, selection)       # ligand itself always opaque

    # Non-molecular objects (gaussian surfaces, maps, CGO, measurements...)
    for obj in cmd.get_names("objects"):
        obj_type = cmd.get_type(obj)
        if obj_type == "object:molecule":
            continue
        if obj in white_objects:
            cmd.color("white", obj)
        elif obj_type in _OUTLINE_RECOLORABLE:
            cmd.color("black", obj)
        else:
            # CGO / measurements carry their own colors: drop them.
            cmd.disable(obj)

    cmd.bg_color("black")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 0)
    cmd.set("ray_shadow", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("depth_cue", 0)
    cmd.set("ambient", 1)
    cmd.set("direct", 0)
    cmd.set("reflect", 0)
    cmd.set("specular", 0)
    cmd.set("spec_direct", 0)
    cmd.set_view(view)
    cmd.png(path, width, height, ray=1)


def _outline_alpha(mask, thickness, offset):
    """Anti-aliased outline band (0-1) from a boolean mask."""
    import numpy as np
    from PIL import Image, ImageChops, ImageFilter

    try:
        from scipy import ndimage as ndi
        sd = (ndi.distance_transform_edt(~mask)
              - ndi.distance_transform_edt(mask))
        return np.clip(thickness / 2 + 0.5 - np.abs(sd - offset), 0, 1)
    except ImportError:
        # PIL fallback: square kernel, no sub-pixel precision.
        m = Image.fromarray(mask.astype(np.uint8) * 255)
        k = 2 * int(thickness) + 1
        if offset < 0:
            e = ImageChops.subtract(m, m.filter(ImageFilter.MinFilter(k)))
        elif offset > 0:
            e = ImageChops.subtract(m.filter(ImageFilter.MaxFilter(k)), m)
        else:
            k = int(thickness) + 1 + (int(thickness) % 2 == 0)
            e = ImageChops.subtract(m.filter(ImageFilter.MaxFilter(k)),
                                    m.filter(ImageFilter.MinFilter(k)))
        return np.array(e) / 255.0


# --- PUBLIC COMMANDS ---
@cmd.extend
def outline_ligand(
    selection="organic",
    filename="outline.png",
    mode="seen",
    width=1600,
    height=1200,
    dpi=300,
    thickness=6,
    offset=0,
    color="black",
    seen_min=0.15,
    glass_fade="proportional",
    extra_white="",
    keep_temp=0,
):
    """
    Render the current scene with a clean outline around one selection.

    The scene is ray-traced twice from the same camera: once normally,
    once as a flat black/white visibility mask of the selection. The
    outline is extracted from the mask and composited on the first
    render. The session is restored afterwards, even on error.

    Parameters
    ----------
    selection : str
        Atoms to outline (e.g. "organic", "resn UQ* and resi 4980").

    filename : str
        Output PNG.

    mode : "visible", "seen" or "full"
        visible : outline only where nothing is in front of the ligand
                  (transparent reps count as opaque).
        seen    : outline what the eye sees, including through
                  transparent reps; hidden behind opaque protein.
        full    : complete silhouette, even behind the protein.

    width, height : int
        Image size in pixels (0, 0 = current viewport).

    dpi : int
        DPI written in the beauty render.

    thickness : float
        Outline thickness in pixels.

    offset : float
        <0 inside the silhouette, >0 outside, 0 centered on the edge.

    color : str
        Outline color (PyMOL color name or 0xRRGGBB).

    seen_min : float
        "seen" only: minimum ligand visibility (0-1) to get an outline.

    glass_fade : "proportional" or float 0-1
        "seen" only: outline opacity behind transparent reps.
        "proportional" = equal to the ligand visibility; 1 = no fade.

    extra_white : str
        Space/comma separated non-molecular objects to treat as part of
        the ligand (e.g. a gaussian surface built from it).

    keep_temp : bool
        Keep the temporary folder (beauty / mask renders) for inspection.

    Examples
    --------
        outline_ligand organic, fig.png

        outline_ligand resn UQ* and resi 4980, fig.png, mode=visible, \
            thickness=4, offset=2, color=cb_vermillion

        outline_ligand organic, fig.png, glass_fade=0.4, \
            extra_white=lig_gaussian

    Notes
    -----
    Requires numpy and Pillow; scipy is optional (smoother, uniform
    outline). CGO and measurement objects are ignored in the mask.
    Gaussian surfaces are isosurface objects, so they're correctly painted black (occluding) in the mask, and their object-level transparency is respected in "seen" mode. If you build a gaussian surface of the ligand, pass it via extra_white=lig_gaussian so it counts as ligand.
    
    CGO (add_axes) and measurement objects (show_bonds) carry their own colors and can't be reliably recolored, so they're disabled in the mask pass — they still appear in the final image, just never occlude the outline. If you want the axes left out of the figure entirely, disable them before calling.
    """
    import numpy as np
    from PIL import Image, ImageFilter

    selection = str(selection).strip()
    if cmd.count_atoms(selection) == 0:
        raise ValueError(f"No atoms found in selection: {selection!r}")

    mode = str(mode).strip().lower()
    if mode not in {"visible", "seen", "full"}:
        raise ValueError("mode must be 'visible', 'seen' or 'full'")

    width, height, dpi = int(width), int(height), int(dpi)
    thickness, offset = float(thickness), float(offset)
    seen_min = float(seen_min)
    keep_temp = str(keep_temp).lower() not in {"0", "false", "off", ""}

    glass_fade = str(glass_fade).strip().lower()
    if glass_fade != "proportional":
        glass_fade = float(glass_fade)

    rgb = tuple(
        int(round(255 * c))
        for c in cmd.get_color_tuple(cmd.get_color_index(str(color)))
    )
    white_objects = {
        o for o in re.split(r"[\s,]+", str(extra_white).strip()) if o
    }

    tmp = tempfile.mkdtemp(prefix="pymol_outline_")
    sess_f = os.path.join(tmp, "backup.pse")
    beauty_f = os.path.join(tmp, "beauty.png")
    mask_f = os.path.join(tmp, "mask.png")
    not_lig = f"not ({selection})"

    cmd.save(sess_f)
    view = cmd.get_view()

    try:
        # 1) Beauty render
        cmd.png(beauty_f, width, height, dpi=dpi, ray=1)

        # 2) Mask render = per-pixel ligand visibility
        if mode == "full":
            cmd.hide("everything", not_lig)
            for obj in cmd.get_names("objects"):
                if (cmd.get_type(obj) != "object:molecule"
                        and obj not in white_objects):
                    cmd.disable(obj)
        _outline_flat_render(
            mask_f, selection, white_objects, width, height,
            keep_transparency=(mode == "seen"), view=view,
        )
        vis = np.array(Image.open(mask_f).convert("RGB")).min(axis=2) / 255.0
        mask = vis > (seen_min if mode == "seen" else 0.8)

        # 3) Outline band, faded behind glass in "seen" mode
        alpha = _outline_alpha(mask, thickness, offset)
        if mode == "seen" and glass_fade != 1:
            k = 2 * int(thickness + abs(offset)) + 3
            vis_img = Image.fromarray((vis * 255).astype(np.uint8))
            vis_band = np.array(vis_img.filter(ImageFilter.MaxFilter(k))) / 255.0
            if glass_fade == "proportional":
                fade = vis_band
            else:
                fade = np.where(vis_band < 0.98, glass_fade, 1.0)
            alpha = alpha * fade

        layer = np.zeros((*mask.shape, 4), dtype=np.uint8)
        layer[..., :3] = rgb
        layer[..., 3] = (alpha * 255).astype(np.uint8)

        # 4) Composite
        beauty = Image.open(beauty_f).convert("RGBA")
        Image.alpha_composite(beauty, Image.fromarray(layer, "RGBA")).save(filename)
        print(f" Outline saved: {os.path.abspath(filename)}")
        if keep_temp:
            print(f" Temp renders kept in: {tmp}")

    finally:
        cmd.load(sess_f)
        if not keep_temp:
            for f in (sess_f, beauty_f, mask_f):
                if os.path.exists(f):
                    os.remove(f)
            os.rmdir(tmp)

    return os.path.abspath(filename)

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
