"""
make_ovito_session.py - builds a ready-made OVITO session (woven_view.ovito).

Run once after (or during) the simulation:
    python make_ovito_session.py              # uses dump.woven
Then on Windows just double-click  D:\\extras\\woven_view.ovito
(or OVITO -> File -> Load Session State). Everything is already set up:
  - warp yarns as solid purple surfaces
  - weft yarns as solid green surfaces
  - broken-bond beads (damage / crack) as red spheres
  - translucent matrix (switch on with the eye icon in the pipeline list)
Just drag the time slider or press Play.

If the dump has per-bead stress (sxx, von_mises columns from in.tensile_sheet),
a second session  woven_stress.ovito  is written: yarns coloured by the axial
stress sigma_xx (blue = compression, red = tension). To show von Mises instead,
pick "von_mises" in the Color coding modifier's property list.

Keep the .ovito files in the same folder as dump.woven.
Needs: pip install ovito   (OVITO desktop must be the same version or newer)
"""
import sys
import ovito
from ovito.io import import_file
import numpy as np
from ovito.modifiers import (SelectTypeModifier, InvertSelectionModifier, ColorCodingModifier,
                             DeleteSelectedModifier, ConstructSurfaceModifier,
                             FreezePropertyModifier, ExpressionSelectionModifier,
                             AssignColorModifier, ComputePropertyModifier)

fname   = sys.argv[1] if len(sys.argv) > 1 else "dump.woven"
session = "woven_view.ovito"
WARP, WEFT, MATRIX, DAMAGE = (0.42, 0.36, 0.85), (0.55, 0.90, 0.65), (0.8, 0.8, 0.82), (0.9, 0.1, 0.1)

def keep_type(p, t):
    p.modifiers.append(SelectTypeModifier(property="Particle Type", types={t}))
    p.modifiers.append(InvertSelectionModifier())
    p.modifiers.append(DeleteSelectedModifier())

def solid(t, color, transparency=0.0, title=""):
    p = import_file(fname)
    keep_type(p, t)
    s = ConstructSurfaceModifier(radius=1.6, smoothing_level=12)
    s.vis.surface_color = color
    s.vis.surface_transparency = transparency
    s.vis.show_cap = False
    p.modifiers.append(s)
    p.source.data.particles.vis.enabled = False   # show surface only, no points
    p.source.data.cell.vis.enabled = False
    p.add_to_scene()
    return p

solid(1, WARP)
solid(2, WEFT)

# damage: beads that have lost bonds since frame 0 (needs c_nb in the dump),
# coloured by HOW MANY bonds they lost (yellow = 1, dark red = many)
probe = import_file(fname)
has_nb = "c_nb" in probe.source.data.particles

def add_lost(p):
    p.modifiers.append(FreezePropertyModifier(source_property="c_nb",
                                              destination_property="nb0", freeze_at=0))
    p.modifiers.append(ComputePropertyModifier(output_property="bonds_lost",
                                               expressions=["nb0 - c_nb"]))

if has_nb:
    d = import_file(fname)
    add_lost(d)
    d.modifiers.append(ExpressionSelectionModifier(expression="bonds_lost < 1"))
    d.modifiers.append(DeleteSelectedModifier())
    d.modifiers.append(ColorCodingModifier(property="bonds_lost", start_value=1, end_value=8,
                                           gradient=ColorCodingModifier.Hot()))
    d.source.data.particles.vis.radius = 0.6
    d.source.data.cell.vis.enabled = False
    d.add_to_scene()
else:
    print("note: no c_nb column in", fname, "- damage pipeline skipped")

m = solid(3, MATRIX, transparency=0.75)
m.modifiers[-1].vis.enabled = False     # matrix off by default (tick it on in OVITO)

ovito.scene.save(session)
nfr = probe.source.num_frames
print(f"wrote {session} with {len(ovito.scene.pipelines)} pipelines "
      f"({nfr} frames). Double-click it in OVITO.")

# ---------------- damage map session: solid yarns coloured by bonds lost -----
if has_nb:
    for p in list(ovito.scene.pipelines):
        p.remove_from_scene()
    for t in (1, 2):
        p = import_file(fname)
        add_lost(p)
        keep_type(p, t)
        s_ = ConstructSurfaceModifier(radius=1.6, smoothing_level=12, transfer_properties=True)
        s_.vis.show_cap = False
        p.modifiers.append(s_)
        p.modifiers.append(ColorCodingModifier(operate_on="surface_vertices", property="bonds_lost",
                                               start_value=0, end_value=4,
                                               gradient=ColorCodingModifier.Jet()))
        p.source.data.particles.vis.enabled = False
        p.source.data.cell.vis.enabled = False
        p.add_to_scene()
    ovito.scene.save("woven_damage.ovito")
    print("wrote woven_damage.ovito (yarns coloured by bonds lost: blue = intact, red = broken)")

# ---------------- stress map session ---------------------------------------
if "sxx" in probe.source.data.particles:
    for p in list(ovito.scene.pipelines):
        p.remove_from_scene()
    # fixed colour range, same for all frames: 98th percentile of |sxx| in yarns
    vals = []
    for fr in sorted(set(np.linspace(0, nfr - 1, min(nfr, 8)).astype(int))):
        dd = probe.compute(fr)
        t = np.asarray(dd.particles["Particle Type"])
        vals.append(np.abs(np.asarray(dd.particles["sxx"])[t < 3]))
    vmax = float(np.percentile(np.concatenate(vals), 98)) or 1.0

    def stress_solid(t, transparency=0.0):
        p = import_file(fname)
        keep_type(p, t)
        s = ConstructSurfaceModifier(radius=1.6, smoothing_level=12, transfer_properties=True)
        s.vis.surface_transparency = transparency
        s.vis.show_cap = False
        p.modifiers.append(s)
        p.modifiers.append(ColorCodingModifier(operate_on="surface_vertices", property="sxx",
                                               start_value=-vmax, end_value=vmax,
                                               gradient=ColorCodingModifier.BlueWhiteRed()))
        p.source.data.particles.vis.enabled = False
        p.source.data.cell.vis.enabled = False
        p.add_to_scene()
        return p

    stress_solid(1)
    stress_solid(2)
    m = stress_solid(3, transparency=0.7)
    m.modifiers[-2].vis.enabled = False      # matrix off by default
    ovito.scene.save("woven_stress.ovito")
    print(f"wrote woven_stress.ovito (sigma_xx colour range +/-{vmax:.3g} LJ)")
else:
    print("note: no sxx column in", fname, "- rerun with the updated in.tensile_sheet for stress maps")
