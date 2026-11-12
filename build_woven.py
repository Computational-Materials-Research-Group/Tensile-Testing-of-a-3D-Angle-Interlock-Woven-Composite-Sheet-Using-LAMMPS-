"""
build_woven.py - LAMMPS data file for a thin 3D angle-interlock woven
composite sheet (warp + weft yarns + matrix), coarse-grained bonded-particle
model in LJ reduced units (bead spacing d = 1 sigma).

  periodic = False : finite rectangular sheet with free edges (tensile coupon)
                     -> run with in.tensile_sheet
  periodic = True  : periodic in-plane RVE (infinite sheet)
                     -> run with in.tensile_periodic

Bead types : 1 = warp yarn, 2 = weft yarn, 3 = matrix
Bond types : 1 = fibre axial, 2 = intra-yarn lateral,
             3 = matrix-matrix, 4 = yarn-matrix interface
"""
import numpy as np
from scipy.spatial import cKDTree

# ======================= INPUTS ===========================================
periodic   = False   # False = finite sheet with free edges
d          = 1.0     # bead spacing (= sigma)
n_layers   = 2       # weft layers through thickness (>= 2; 2 = thinnest sheet)
n_rep_x    = 8       # weave repeats along x (length, loading direction)
n_warp_y   = 8       # warp yarn columns along y (width)
px         = 14.0    # weft spacing along x (= half warp wavelength)
py         = 10.0    # warp spacing along y
weft_a, weft_b = 4.0, 1.5   # weft ellipse semi-axes: width (x), thickness (z)
warp_a, warp_b = 4.0, 1.5   # warp ellipse semi-axes: width (y), thickness
gap        = 0.5     # clearance between yarns
beta       = 2.5     # crimp shape: 0 = sinusoid, larger = flat over wefts
r_bond     = 1.5     # bond search radius (face + diagonal neighbours)
r_excl     = 0.9     # matrix beads closer than this to a yarn bead are removed
add_matrix = True
outfile    = "woven.data"
# ==========================================================================

assert n_layers >= 2, "need at least 2 weft layers for the warp to interlock"
Lx, Ly = 2 * n_rep_x * px, n_warp_y * py
h  = weft_b + warp_b + 0.5 * gap          # warp-centre offset above/below a weft
dz = 2 * h + gap                           # weft layer spacing
A  = 0.5 * dz + h                          # warp undulation amplitude
x0 = 0.0 if periodic else 0.5 * px         # finite sheet: keep wefts inside box

def ellipse_hex(a, b, s):
    """hexagonally packed filament positions inside an ellipse"""
    pts, rs = [], s * np.sqrt(3) / 2
    for iv in range(-int(b / rs) - 1, int(b / rs) + 2):
        v = iv * rs
        for iu in range(-int(a / s) - 2, int(a / s) + 3):
            u = iu * s + 0.5 * s * (iv % 2)
            if (u / a) ** 2 + (v / b) ** 2 <= 1.0:
                pts.append((u, v))
    return np.array(pts)

pos, typ, mol, fil, idx, nfil = [], [], [], [], [], {}
yarn_id = fil_id = 0

# ---- weft yarns: straight along y -----------------------------------------
sec = ellipse_hex(weft_a, weft_b, d)
nb  = int(round(Ly / d))
for k in range(n_layers):
    for i in range(2 * n_rep_x):
        yarn_id += 1
        for (u, v) in sec:
            fil_id += 1; nfil[fil_id] = nb
            for n in range(nb):
                y = n * Ly / nb if periodic else (n + 0.5) * Ly / nb
                pos.append((x0 + i * px + u, y, k * dz + v))
                typ.append(2); mol.append(yarn_id); fil.append(fil_id); idx.append(n)

# ---- warp yarns: undulate along x, interlocking adjacent weft layers -------
sec = ellipse_hex(warp_a, warp_b, d)
xs  = np.linspace(0.0, Lx, 8001)
for j in range(n_warp_y):
    yj = (j + 0.5) * py
    if n_layers == 2:                      # one channel, phase alternates per column
        channels, phase = [0], np.pi * (j % 2)
    else:                                  # alternate channels per column
        channels, phase = range(j % 2, n_layers - 1, 2), 0.0
    c  = np.cos(np.pi * (xs - x0) / px + phase)   # crests aligned with wefts
    sn = np.sin(np.pi * (xs - x0) / px + phase)
    if beta > 0:
        shape = A * np.tanh(beta * c) / np.tanh(beta)
        dzdx  = -A * beta * np.pi / px * sn / np.cosh(beta * c) ** 2 / np.tanh(beta)
    else:
        shape, dzdx = A * c, -A * np.pi / px * sn
    nrm = np.sqrt(1 + dzdx ** 2)
    nx_, nz_ = -dzdx / nrm, 1.0 / nrm
    for m in channels:
        yarn_id += 1
        zc = (m + 0.5) * dz + shape
        for (u, v) in sec:
            fil_id += 1
            X, Z = xs + v * nx_, zc + v * nz_
            s  = np.concatenate([[0], np.cumsum(np.hypot(np.diff(X), np.diff(Z)))])
            nb = int(round(s[-1] / d)); nfil[fil_id] = nb
            sk = np.arange(nb) * s[-1] / nb if periodic else (np.arange(nb) + 0.5) * s[-1] / nb
            for n, (xk, zk) in enumerate(zip(np.interp(sk, s, X), np.interp(sk, s, Z))):
                pos.append((xk, yj + u, zk))
                typ.append(1); mol.append(yarn_id); fil.append(fil_id); idx.append(n)

pos = np.array(pos); typ = np.array(typ); mol = np.array(mol)
fil = np.array(fil); idx = np.array(idx)
if periodic:
    pos[:, 0] %= Lx; pos[:, 1] %= Ly
    pos[pos[:, 0] >= Lx, 0] = 0.0; pos[pos[:, 1] >= Ly, 1] = 0.0
nyarn = len(pos)
zlo, zhi = pos[:, 2].min() - 1.0, pos[:, 2].max() + 1.0

def tree(p):
    if periodic:
        q = p.copy(); q[:, 2] -= zlo
        return cKDTree(q, boxsize=[Lx, Ly, 3 * (zhi - zlo) + 10.0])
    return cKDTree(p)
def shift(p):
    return np.column_stack([p[:, 0], p[:, 1], p[:, 2] - zlo]) if periodic else p

ty = tree(pos)
close = ty.query_pairs(0.8 * d, output_type="ndarray")
cross = close[mol[close[:, 0]] != mol[close[:, 1]]]
print(f"yarn beads: {nyarn}, yarns: {yarn_id}, inter-yarn overlaps (<0.8d): {len(cross)}")
if len(cross):
    print("  WARNING: increase 'gap', 'px' or 'beta', or reduce yarn thickness")

# ---- matrix: simple-cubic lattice filling the space between yarns ----------
if add_matrix:
    nxm, nym = int(round(Lx / d)), int(round(Ly / d))
    nzm = int((zhi - zlo - 1) / d)
    g = np.stack(np.meshgrid((np.arange(nxm) + 0.5) * Lx / nxm,
                             (np.arange(nym) + 0.5) * Ly / nym,
                             zlo + 0.5 + np.arange(nzm) * d, indexing="ij"), -1).reshape(-1, 3)
    dist, _ = ty.query(shift(g))
    g = g[dist > r_excl]
    pos = np.vstack([pos, g]); typ = np.concatenate([typ, np.full(len(g), 3)])
    mol = np.concatenate([mol, np.zeros(len(g), int)])
    fil = np.concatenate([fil, np.zeros(len(g), int)])
    idx = np.concatenate([idx, np.zeros(len(g), int)])
print(f"total beads: {len(pos)}  (matrix: {len(pos) - nyarn})")

# ---- bonds -----------------------------------------------------------------
pairs = tree(pos).query_pairs(r_bond, output_type="ndarray")
a, b = pairs[:, 0], pairs[:, 1]
ya, yb = typ[a] < 3, typ[b] < 3
nb_f = np.array([nfil.get(f, 1) for f in fil[a]])
dn = np.abs(idx[a] - idx[b])
neighbour = (dn == 1) | ((dn == nb_f - 1) & periodic)
axial = ya & (fil[a] == fil[b]) & neighbour
btype = np.zeros(len(pairs), int)
btype[axial] = 1
btype[~axial & ya & yb & (mol[a] == mol[b])] = 2
btype[~ya & ~yb] = 3
btype[ya ^ yb] = 4
keep = btype > 0                 # different yarns: contact only, no bond
a, b, btype = a[keep], b[keep], btype[keep]
for t in range(1, 5):
    print(f"bond type {t}: {np.sum(btype == t)}")

# ---- write -----------------------------------------------------------------
if periodic:
    xb, yb_ = (0.0, Lx), (0.0, Ly)
else:
    xb, yb_ = (pos[:, 0].min() - 1, pos[:, 0].max() + 1), (pos[:, 1].min() - 1, pos[:, 1].max() + 1)
with open(outfile, "w") as f:
    f.write("3D angle-interlock woven composite (bonded-particle model)\n\n")
    f.write(f"{len(pos)} atoms\n{len(a)} bonds\n3 atom types\n4 bond types\n\n")
    f.write(f"{xb[0]:.6f} {xb[1]:.6f} xlo xhi\n{yb_[0]:.6f} {yb_[1]:.6f} ylo yhi\n")
    f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")
    f.write("Masses\n\n1 1.0\n2 1.0\n3 1.0\n\nAtoms # bond\n\n")
    f.writelines(f"{i+1} {mo} {t} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n"
                 for i, (p, t, mo) in enumerate(zip(pos, typ, mol)))
    f.write("\nBonds\n\n")
    f.writelines(f"{i+1} {t} {p+1} {q+1}\n" for i, (t, p, q) in enumerate(zip(btype, a, b)))
print(f"wrote {outfile}: {'periodic RVE' if periodic else 'finite sheet'} "
      f"{xb[1]-xb[0]:.1f} x {yb_[1]-yb_[0]:.1f} x {zhi-zlo:.1f} (L x W x t)")
