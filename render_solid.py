"""
render_solid.py - render the woven composite as solid yarns (not points) with OVITO.

Usage:
  python render_solid.py                      # woven.data  -> woven_solid.png
  python render_solid.py dump.woven           # last frame of the tensile run
  python render_solid.py dump.woven 10        # frame 10
  python render_solid.py dump.woven all       # every frame -> movie frames

Needs:  pip install ovito
"""
import sys
from ovito.io import import_file
from ovito.modifiers import (SelectTypeModifier, InvertSelectionModifier,
                             DeleteSelectedModifier, ConstructSurfaceModifier,
                             SliceModifier)
from ovito.vis import Viewport, TachyonRenderer

# ---------------- settings ------------------------------------------------
fname        = sys.argv[1] if len(sys.argv) > 1 else "woven.data"
frame_arg    = sys.argv[2] if len(sys.argv) > 2 else "last"
show_matrix  = False       # True = translucent matrix around the yarns
cut_matrix   = True        # cut the matrix open (half removed) to show the weave
WARP_COLOR   = (0.42, 0.36, 0.85)   # purple
WEFT_COLOR   = (0.55, 0.90, 0.65)   # green
MATRIX_COLOR = (0.80, 0.80, 0.82)
alpha_radius = 1.6         # surface probe radius (bead spacing = 1)
smoothing    = 12          # higher = smoother yarns
size         = (1800, 1000)
# ---------------------------------------------------------------------------

def load():
    if fname.endswith(".data"):
        return import_file(fname, atom_style="bond")
    return import_file(fname)

def hide_particles(frame, data):
    data.particles.vis.enabled = False
    if data.particles.bonds is not None:
        data.particles.bonds.vis.enabled = False
    data.cell.vis.enabled = False        # no box outline

def solid(ptype, color, transparency=0.0, slice_matrix=False):
    p = load()
    p.modifiers.append(SelectTypeModifier(property="Particle Type", types={ptype}))
    p.modifiers.append(InvertSelectionModifier())
    p.modifiers.append(DeleteSelectedModifier())
    if slice_matrix:                     # remove the half of the matrix facing the camera
        p.modifiers.append(SliceModifier(normal=(0, 1, 0), distance=0.0,
                                         slab_width=0.0, inverse=True))
    surf = ConstructSurfaceModifier(radius=alpha_radius, smoothing_level=smoothing)
    surf.vis.surface_color = color
    surf.vis.surface_transparency = transparency
    surf.vis.show_cap = False
    p.modifiers.append(surf)
    p.modifiers.append(hide_particles)
    p.add_to_scene()
    return p

warp = solid(1, WARP_COLOR)
weft = solid(2, WEFT_COLOR)
if show_matrix:
    m = solid(3, MATRIX_COLOR, transparency=0.55, slice_matrix=cut_matrix)
    if cut_matrix:   # place the cut at mid-width of the sheet
        cell = m.source.data.cell
        ymid = cell[1, 3] + 0.5 * cell[1, 1]
        m.modifiers[3].distance = ymid

vp = Viewport(type=Viewport.Type.Perspective, camera_dir=(1.0, 1.3, -0.9), fov=0.35)
vp.zoom_all(size)
renderer = TachyonRenderer(ambient_occlusion=True, shadows=True, antialiasing_samples=8)

nframes = warp.source.num_frames
base = fname.rsplit(".", 1)[0]
if frame_arg == "all":
    vp.render_anim(filename=f"{base}_solid.png", size=size, renderer=renderer,
                   background=(1, 1, 1), range=(0, nframes - 1))
    print(f"wrote {nframes} frames {base}_solid*.png")
else:
    fr = nframes - 1 if frame_arg == "last" else int(frame_arg)
    vp.render_image(filename=f"{base}_solid.png", size=size, frame=fr,
                    renderer=renderer, background=(1, 1, 1))
    print(f"wrote {base}_solid.png (frame {fr} of {nframes})")
