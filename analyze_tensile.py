"""
analyze_tensile.py - results of the woven-sheet tensile test.

Reads stress_strain.txt (written by in.tensile_sheet / in.tensile_periodic)
and reports Young's modulus, tensile strength, strain at peak, failure strain
and bond damage; saves stress_strain.png.

Usage:  python analyze_tensile.py [stress_strain.txt]
Optional real units: set SIGMA_TO_MPA (stress unit in MPa) below.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fname        = sys.argv[1] if len(sys.argv) > 1 else "stress_strain.txt"
fit_limit    = 0.005     # fit modulus on strain <= 0.5 % ...
fit_frac     = 0.4       # ... and stress <= 40 % of peak (linear part only)
SIGMA_TO_MPA = None      # e.g. 50.0 if 1 LJ stress unit = 50 MPa after calibration

d = np.loadtxt(fname, comments="#")
strain, stress = d[:, 1], d[:, 2]
bonds = d[:, -1]
unit = "LJ"
if SIGMA_TO_MPA:
    stress = stress * SIGMA_TO_MPA; unit = "MPa"

ipk = np.argmax(stress)
uts, e_uts = stress[ipk], strain[ipk]
lin = (strain > 1e-5) & (strain <= fit_limit) & (stress <= fit_frac * uts) & (np.arange(len(strain)) < ipk)
if lin.sum() < 3:                                   # fall back to first 10 % of points before peak
    lin = np.arange(len(strain)) < max(3, ipk // 10)
E, c0 = np.polyfit(strain[lin], stress[lin], 1)

after = np.where((np.arange(len(stress)) > ipk) & (stress < 0.1 * uts))[0]
e_fail = strain[after[0]] if len(after) else None
broken = 1.0 - bonds / bonds[0]
first_break = strain[np.argmax(bonds < bonds[0])] if np.any(bonds < bonds[0]) else None
toughness = np.trapezoid(np.clip(stress, 0, None), strain) if hasattr(np, "trapezoid") \
            else np.trapz(np.clip(stress, 0, None), strain)

print(f"data points            : {len(strain)}  (max strain {strain.max():.4f})")
print(f"Young's modulus E      : {E:.4g} {unit}")
print(f"tensile strength (UTS) : {uts:.4g} {unit}")
print(f"strain at UTS          : {e_uts:.4f}  ({100*e_uts:.2f} %)")
print(f"first bond breaks at   : {first_break if first_break is None else f'{first_break:.4f}'}")
print(f"failure strain (<10%UTS): {'not reached yet' if e_fail is None else f'{e_fail:.4f}'}")
print(f"bonds broken at end    : {100*broken[-1]:.2f} %")
print(f"toughness (area)       : {toughness:.4g} {unit}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(100 * strain, stress, lw=1.5, color="#5b4fd6")
xs = np.linspace(0, strain[lin].max() * 1.6, 20)
ax[0].plot(100 * xs, E * xs + c0, "--", color="gray", label=f"E = {E:.3g} {unit}")
ax[0].plot(100 * e_uts, uts, "o", color="crimson", label=f"UTS = {uts:.3g} {unit}")
ax[0].set_xlabel("engineering strain (%)"); ax[0].set_ylabel(f"stress ({unit})")
ax[0].set_title("Stress–strain"); ax[0].legend(); ax[0].grid(alpha=0.3)
ax[1].plot(100 * strain, 100 * broken, color="#3aa86b", lw=1.5)
ax[1].set_xlabel("engineering strain (%)"); ax[1].set_ylabel("bonds broken (%)")
ax[1].set_title("Damage evolution"); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.savefig("stress_strain.png", dpi=200)
print("saved stress_strain.png")
