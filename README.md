# Tensile Testing of a 3D Angle-Interlock Woven Composite Sheet Using LAMMPS (BPM Bonded-Particle Model)


![LAMMPS](https://img.shields.io/badge/LAMMPS-BPM%20Package-blue?style=for-the-badge)
![Model](https://img.shields.io/badge/Model-3D%20Woven%20Composite-purple?style=for-the-badge)
![Bonds](https://img.shields.io/badge/bond__style-bpm%2Fspring-orange?style=for-the-badge)
![Loading](https://img.shields.io/badge/Loading-Displacement%20Controlled-black?style=for-the-badge)
![OVITO](https://img.shields.io/badge/OVITO-Visualization-green?style=for-the-badge)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24155083-blue?style=for-the-badge)](https://doi.org/10.5281/zenodo.24155083)

A **coarse-grained bonded-particle simulation of the uniaxial tensile test of a thin 3D angle-interlock woven composite sheet**, using LAMMPS and the BPM package. Warp and weft yarns are bundles of bead filaments with elliptical cross-sections, the warp yarns undulate through the thickness and lock two layers of weft yarns together, and a bead matrix fills the space between the yarns. Every bond is a breakable spring (`bond_style bpm/spring`). Both ends of the sheet are clamped in rigid grips; the left grip is fixed and the right grip is pulled at constant speed until the sheet fails.

The simulation captures **elastic loading of the crimped woven architecture**, **straightening of the warp crimp**, **the onset and growth of fibre, matrix and interface damage**, and **the loss of load-carrying capacity**. Strain, stress, grip force and the number of intact bonds are written every 200 steps, and the full trajectory (with damage and per-bead stress) is written for visualization in OVITO, where the yarns are rendered as solid bodies.


<img width="640" height="480" alt="composite woven" src="https://github.com/user-attachments/assets/a831263c-6c9a-4a1d-a550-430ae4ab9977" />


---

## Physics

The simulation models the composite as a network of breakable springs between beads:

- Each yarn is a bundle of 21 filaments, hexagonally packed in an elliptical cross-section; each filament is a chain of beads one bead spacing apart
- Weft yarns are straight (along y) and stacked in two layers; warp yarns (along x, the loading direction) undulate between the layers and interlock them
- The warp path is a smoothed square wave: flat over a weft, steep between wefts, like real crimp; adjacent warp columns are half a wavelength out of phase
- The matrix is a simple-cubic bead lattice that fills the space between the yarns
- Bonds join every pair of beads within 1.5 bead spacings (face and diagonal neighbours), except between different yarns: yarns touch each other only by contact and are held together by the matrix
- Four bond types: fibre axial, fibre transverse (within a yarn), matrix, and yarn/matrix interface, each with its own stiffness and failure strain
- Bonds are springs that resist stretching only; bending stiffness of a yarn comes from the bundle of cross-linked filaments, not from beam bonds
- Each bond stores its initial length as the reference, so the structure starts stress-free
- A bond breaks when its strain exceeds its critical strain; broken bonds become contact pairs (`pair_style bpm/spring`)
- Energy is dissipated by velocity damping inside the bonds and contacts, which does not resist rigid motion or uniform stretching

---

## Loading Regime

| Quantity                       | Value                    |
| ------------------------------ | ------------------------ |
| Nominal strain rate (`erate`)  | 5.0e-5 per time          |
| Gauge length                   | 205.42                   |
| Grip speed (`vpull`)           | 1.03e-2                  |
| Timestep                       | 0.005                    |
| Loading steps                  | 400,000                  |
| Final engineering strain       | 10 %                     |
| Strain per dump frame          | 0.1 %                    |
| Thermostat during loading      | none (bond damping only) |

All quantities are in non-dimensional (`lj`) units. At the start of loading the gauge section is given a linear velocity profile (0 at the left grip, `vpull` at the right grip), so the test starts without a stress wave from the grip.

---

## Model and Bond Parameters

| Parameter                     | Value          | Role                                                     |
| ----------------------------- | -------------- | -------------------------------------------------------- |
| Bead spacing `d`              | 1.0            | Length unit of the model                                 |
| Weft layers `n_layers`        | 2              | Through-thickness layers (thinnest interlock)            |
| Weave repeats `n_rep_x`       | 8              | Repeats along the loading direction (length)             |
| Warp columns `n_warp_y`       | 8              | Warp yarns across the width                              |
| Weft spacing `px`             | 14.0           | Half the warp wavelength                                 |
| Warp spacing `py`             | 10.0           | Spacing of warp columns                                  |
| Weft semi-axes `a`, `b`       | 4.0, 1.5       | Width (x) and thickness (z) of a weft                    |
| Warp semi-axes `a`, `b`       | 4.0, 1.5       | Width (y) and thickness of a warp                        |
| Yarn clearance `gap`          | 0.5            | Space between neighbouring yarns                         |
| Crimp shape `beta`            | 2.5            | 0 = sinusoid, larger = flatter over the wefts            |
| Bond search radius            | 1.5            | Face and diagonal neighbours are bonded                  |
| Matrix exclusion radius       | 0.9            | Matrix beads closer than this to a yarn bead are removed |
| Bead mass                     | 1.0            | All three bead types                                     |
| Grip length                   | 10.0           | Clamped length at each end                               |

### Bond Coefficients (`bond_style bpm/spring`)

| Bond type | Meaning                 | Stiffness `k` | Critical strain `eps_c` | Damping `gamma` |
| --------- | ----------------------- | ------------- | ----------------------- | --------------- |
| 1         | Fibre, axial            | 100.0         | 0.030                   | 1.0             |
| 2         | Fibre, transverse       | 20.0          | 0.050                   | 1.0             |
| 3         | Matrix                  | 5.0           | 0.080                   | 1.0             |
| 4         | Yarn/matrix interface   | 3.0           | 0.040                   | 1.0             |

Contact between non-bonded beads: `pair_coeff * * 1.0 1.0 1.0` (stiffness 1.0, contact distance 1.0, damping 1.0).

The geometry is generated by `build_woven.py` and stored in `woven.data`.

---

## Geometry

```
Side view through one warp yarn (x-z plane, units: lj, non-dimensional):

            warp (purple) over weft     warp under weft
                 ________                                  ________
z = +11.8      /        \                                /        \
              /  ( W1 )  \           ( W1 )             /  ( W1 )  \      <- weft layer 1, z = 7
             |            |                            |            |
             |   ( W0 )   |          ( W0 )            |   ( W0 )   |     <- weft layer 0, z = 0
z = -4.8                   \________________________/
         |<---- px = 14 ---->|

  the warp lies on top of layer 1 at one weft column and under layer 0 at the next:
  it passes through both layers and locks them together (layer-to-layer angle interlock)

Top view: 8 warp columns (y spacing 10) cross 16 weft columns (x spacing 14);
neighbouring warp columns are half a wavelength out of phase.

Tensile coupon (x = loading direction):

  +--------+------------------------------------------------+--------+
  |  LEFT  |                 gauge length 205.4             | RIGHT  |
  |  GRIP  |                                                |  GRIP  |  ---> v
  | fixed  |                                                | pulled |
  +--------+------------------------------------------------+--------+
   10.0                                                       10.0

x, y, z: all non-periodic (boundary s s s), free edges and surfaces
```

- **Sheet (data box)**: 227.4 x 81.0 x 17.2 (length x width x thickness)
- **Weave**: 2 weft layers 7.0 apart, warp undulation amplitude 6.75, warp offset 3.25 above or below a weft centre
- **Warp path**: smallest radius of curvature 2.67 (larger than the half-thickness 1.5 of the warp, so filaments do not fold), steepest slope 75 degrees
- **Overlap check**: 0 bead pairs from different yarns closer than 0.8 (reported by `build_woven.py`)
- **Measured at the start of loading**: gauge length 205.42, width 79.0, thickness 15.73, cross-section 1242.82
- **Yarn bead fraction**: 45.3 % of all beads

### System Size

| Component                       | Count       |
| ------------------------------- | ----------- |
| Beads in total                  | 251,408     |
| of which warp (type 1)          | 60,144      |
| of which weft (type 2)          | 53,760      |
| of which matrix (type 3)        | 137,504     |
| Yarns (8 warp + 32 weft)        | 40          |
| Filaments per yarn              | 21          |
| Bonds in total                  | 1,764,496   |

---

## Atom and Bond Types

| Type | Kind | Role                                                    |
| ---- | ---- | ------------------------------------------------------- |
| 1    | Atom | Warp yarn bead (along x, the loading direction)         |
| 2    | Atom | Weft yarn bead (along y)                                |
| 3    | Atom | Matrix bead                                             |

| Bond type | Kind                     | Count     |
| --------- | ------------------------ | --------- |
| 1         | Fibre axial              | 113,064   |
| 2         | Fibre transverse         | 635,520   |
| 3         | Matrix-matrix            | 916,168   |
| 4         | Yarn-matrix interface    | 99,744    |

The molecule ID of a yarn bead is the number of its yarn (1 to 40); matrix beads have molecule ID 0.

---

## Simulation Phases

| Phase     | Run Steps   | Simulated Time | Description                                                          |
| --------- | ----------- | -------------- | -------------------------------------------------------------------- |
| Relax     | 2,000       | 10             | Grips held, light Langevin damping (T = 1e-6), structure at rest     |
| Loading   | 400,000     | 2,000          | Left grip fixed, right grip pulled at constant speed, no thermostat  |
| **Total** | **402,000** | **2,010**      | Timestep 0.005                                                       |

A full run took **40 min 07 s on 14 MPI ranks** (about 146 timesteps per second), with 0 dangerous neighbour-list builds.

---

## Expected Event Timeline

The events below were read from the thermo output of one run with the settings above.

| Event                                           | Strain   | Stress (lj) | Intact bonds | Approx. dump frame |
| ----------------------------------------------- | -------- | ----------- | ------------ | ------------------ |
| Start of loading                                | 0        | 0           | 1,764,496    | 0                  |
| Linear elastic response                         | 0.05 %   | 0.0033      | 1,764,496    | 0 to 1             |
| Linear elastic response                         | 0.10 %   | 0.0066      | 1,764,496    | 1                  |
| Bond failure has started (328 bonds lost)       | 0.90 %   | 0.0578      | 1,764,168    | 9                  |
| Damage growing, curve flattening (2,049 lost)   | 1.05 %   | 0.0625      | 1,762,447    | 10 to 11           |

The initial modulus is **E = 6.6** (lj units). Damage first appears at the yarn crossovers, where the crimped warp straightens and loads the interface. Run `analyze_tensile.py` on your own `stress_strain.txt` for the peak stress, the strain at the peak, the failure strain and the final damage.

---

## Governing Equations

### Weave Geometry

```
h  = b_weft + b_warp + gap/2              offset of a warp centre from a weft centre    = 3.25
dz = 2 h + gap                            spacing of the weft layers                    = 7.0
A  = dz/2 + h                             warp undulation amplitude                     = 6.75

weft layer k      : z = k dz,  x = x0 + i px                 (straight along y)
warp column j     : y = (j + 1/2) py
warp centre line  : z(x) = dz/2 + A tanh(beta cos(pi (x - x0)/px + phi_j)) / tanh(beta)
phase             : phi_j = pi (j mod 2)                     (two weft layers)
filament (u, v)   : r = c(x) + u e_y + v n(x)                n = in-plane normal of the warp path
```

### Bond Force (LAMMPS `bpm/spring`)

```
F = k (r - r0) / r0  -  gamma (n . dv) n           r0 = length when the bond is first computed

the bond breaks when |r - r0| / r0 > eps_c (it is then removed and the pair interacts by contact only)
```

### Loading Kinematics

```
L_g      = 205.42                              gauge length (between the grips)
v_pull   = erate * L_g = 5e-5 * 205.42 = 1.03e-2
strain   = (x_cm(right grip) - x_cm(left grip) - initial distance) / L_g
```

### Measured Stress

```
F      = (F_left - F_right) / 2              sum of x-forces on each grip, averaged
stress = F / A0                              A0 = width x thickness at the start of loading = 1242.82
```

---

## Group Definitions

```
gripL    <- beads with x < xmin + 10            (left grip, held by fix move with zero speed)
gripR    <- beads with x > xmax - 10            (right grip, moved by fix move at v_pull)
grips    <- gripL + gripR
mobile   <- all - grips                          (integrated with fix nve)
```

---

## Output Files

| File                 | Frequency          | Content                                                              |
| -------------------- | ------------------ | -------------------------------------------------------------------- |
| `stress_strain.txt`  | every 200 steps    | 5 columns (see the table below), averages of 20 samples              |
| `dump.woven`         | every 4,000 steps  | id, mol, type, x, y, z, c_nb, sxx, von_mises (about 100 frames)      |
| `log.lammps`         | every 2,000 steps  | thermo output                                                        |

### Columns of `stress_strain.txt`

| Column | Quantity                                       |
| ------ | ---------------------------------------------- |
| 1      | Timestep                                       |
| 2      | Engineering strain (`v_strain`)                |
| 3      | Engineering stress (`v_stress`)                |
| 4      | Grip force (`v_force`)                         |
| 5      | Number of intact bonds (`v_bonds`)             |

### Per-Bead Columns of `dump.woven`

| Column       | Quantity                                                            |
| ------------ | ------------------------------------------------------------------- |
| `c_nb`       | Number of intact bonds of the bead (drops when bonds break)         |
| `sxx`        | Stress along x (times the bead volume, about 1), positive = tension |
| `von_mises`  | Von Mises equivalent stress                                         |

`sxx` and `von_mises` are time averages over the 400 steps before each frame.

### Thermo Output

| Column       | Quantity                                 |
| ------------ | ---------------------------------------- |
| `step`       | Timestep                                 |
| `v_strain`   | Engineering strain                       |
| `v_stress`   | Engineering stress                       |
| `v_force`    | Grip force                               |
| `v_bonds`    | Number of intact bonds                   |
| `ke`         | Kinetic energy                           |

---

## Repository Structure

```
woven_tensile_lammps/
|
|-- build_woven.py          # generates woven.data (sheet or periodic RVE)
|-- in.tensile_sheet        # LAMMPS input: tensile test of the finite sheet with grips (this run)
|-- in.tensile_periodic     # LAMMPS input: tensile test of the periodic RVE (infinite sheet)
|-- README.md               # this file
|
|-- optional helper files/
|   |-- analyze_tensile.py      # modulus, strength, failure strain, damage; saves stress_strain.png
|   |-- make_ovito_session.py   # writes ready-made OVITO sessions (solid yarns, damage, stress)
|   |-- render_solid.py         # renders solid-yarn images or animation frames without the GUI
|
|-- generated on run/
    |-- woven.data
    |-- stress_strain.txt
    |-- dump.woven
    |-- log.lammps
    |-- stress_strain.png
    |-- woven_view.ovito
    |-- woven_damage.ovito
    |-- woven_stress.ovito
```

---

## How to Run

### Requirements

- LAMMPS with the **BPM** and **MOLECULE** packages: <https://www.lammps.org> (the case was run with the 22 Jul 2025 release)
- An MPI runtime that matches the one LAMMPS was built with, for parallel runs
- Python with numpy, scipy and matplotlib
- OVITO desktop and the `ovito` Python module (same version or newer for the desktop app) for visualization: <https://www.ovito.org>

Check that LAMMPS has the BPM package:

```
lmp -h | grep bpm/spring
```

LAMMPS with BPM can be installed with conda:

```
conda install -c conda-forge "lammps=*=*mpi_mpich*" -y
```

### Step 1: Build the sheet

Edit the inputs at the top of `build_woven.py` if needed (`periodic = False` for the sheet), then:

```
python build_woven.py
```

The last line should read `wrote woven.data: finite sheet 227.4 x 81.0 x 17.2 (L x W x t)`, and the overlap check should report 0.

### Step 2: Run the tensile test

With MPI on Ubuntu or WSL, with the files on `D:`:

```
cd /mnt/d/extras
nohup mpirun -np 14 lmp -in in.tensile_sheet > run.out 2>&1 &
```

Use `-np` up to the number of physical cores (`nproc`). Follow the run with `tail -f log.lammps`; the run is finished when `jobs` lists nothing and the log ends with `Total wall time`.

### Step 3: Analyse the stress-strain curve

```
python analyze_tensile.py
```

This prints the Young's modulus, the tensile strength, the strain at the peak, the first bond failure, the failure strain, the broken-bond fraction and the toughness, and saves `stress_strain.png`. Set `SIGMA_TO_MPA` in the script to report stresses in MPa after calibration.

### Step 4: Visualize in OVITO

```
python make_ovito_session.py
```

This writes up to three session files next to `dump.woven`. Open one in OVITO with `File > Load Session State` (not `Load File`), then press Play:

| Session file          | What it shows                                                                   |
| --------------------- | ------------------------------------------------------------------------------- |
| `woven_view.ovito`    | Solid warp (purple) and weft (green) yarns; damaged beads coloured by bonds lost |
| `woven_damage.ovito`  | Solid yarns coloured by damage: blue = intact, red = broken                     |
| `woven_stress.ovito`  | Solid yarns coloured by `sxx`: blue = compression, white = 0, red = tension     |

To render images without the GUI:

```
python render_solid.py dump.woven last     # final frame -> dump_solid.png
python render_solid.py dump.woven all      # all frames, for an animation
ffmpeg -framerate 10 -i dump_solid%04d.png -pix_fmt yuv420p tensile.mp4
```

---

## Important: Damping and the Start of Loading

Two settings strongly affect the measured stress, and both were found while testing this case:

| Setting                                                  | Effect                                                                 |
| -------------------------------------------------------- | ---------------------------------------------------------------------- |
| Langevin thermostat kept on during loading               | Drags every bead against its own motion; false bond failure near the moving grip at 0.25 % strain and a flat stress plateau |
| Right grip started suddenly, without a velocity profile   | A stress wave from the grip breaks bonds early                          |
| Langevin off during loading and linear velocity profile (used here) | Linear elastic response, no bonds lost before damage starts |

The Langevin thermostat damps absolute velocity, so in a stretched body it acts like a viscous drag that grows with the distance from the fixed grip. It is used here only for the short relaxation and is switched off (`unfix damp`) before loading; the damping inside the bonds and contacts dissipates energy without resisting uniform stretching. If you change the strain rate, check that the stress does not depend on it by halving `erate`.

---

## What to Look for in Results

### Elastic Loading

The stress rises linearly with the strain and the number of intact bonds stays at 1,764,496. The initial modulus is about 6.6 in lj units.

### Crimp Straightening

The warp yarns carry the load along x. Because they are crimped, they first straighten, which loads the yarn crossovers. In `woven_stress.ovito` the highest tensile stress sits on the crests of the warp yarns where they pass over the wefts.

### Onset of Failure

The first bonds to fail are at the crossovers (interface and transverse fibre bonds, which have the lowest critical strains relative to their local strain). In OVITO the damage appears as vertical bands across the width at the crossover columns.

### Damage Growth and Load Drop

Column 5 of `stress_strain.txt` falls faster as fibres begin to break. After the peak the stress drops while the grip keeps moving; a very low final stress means the warp yarns have broken across the whole width.

### Kinetic Energy

`ke` in the thermo output must stay small compared with the stored elastic energy. A large kinetic energy means the loading is too fast.

---

## Common Errors and Fixes

| Error or symptom                                          | Cause                                                          | Fix                                                                         |
| --------------------------------------------------------- | -------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `D:: command not found`, `cd: D:extras: No such file`     | Windows path typed in a Linux (WSL) terminal                   | Use `cd /mnt/d/extras`                                                      |
| `python: can't open file 'build_woven.py'`                | Command run outside the folder with the files                  | `cd` into the folder that holds the scripts                                 |
| `Unrecognized fix style 'update/special/bonds'`           | Fix removed from recent LAMMPS versions                        | Not needed; it is not used in these inputs                                  |
| Unknown bond or pair style `bpm/spring`                   | LAMMPS built without the BPM package, or too old               | Install a LAMMPS that includes BPM (check with `lmp -h`)                    |
| `WARNING: Inconsistent image flags` (periodic RVE)        | Yarns are continuous across the periodic boundaries            | Harmless                                                                    |
| Overlap warning from `build_woven.py`                     | Yarns too thick or too close for the crimp                     | Increase `gap`, `px` or `beta`, or reduce the yarn thickness                |
| Early bond failure near the moving grip, flat stress      | Langevin thermostat on during loading, or sudden grip start    | Keep `unfix damp` and the velocity ramp of `in.tensile_sheet`               |
| `ModuleNotFoundError: No module named 'ovito'`            | OVITO Python module not installed                              | `pip install ovito`, or the conda package from `conda.ovito.org`            |
| Warning about OVITO from PyPI in a conda environment      | pip OVITO package inside conda                                 | Harmless if the scripts run; otherwise use the conda package                |
| Qt, `libOpenGL` or Vulkan error when rendering            | No graphics libraries on WSL or a server                       | `sudo apt install libopengl0 libegl1 mesa-vulkan-drivers`                   |
| Only points in OVITO                                      | `dump.woven` opened directly                                   | Open a session with `File > Load Session State`                             |
| Damaged beads appear as a cloud of points                 | Many beads have lost bonds after failure                       | Hide the damage pipeline, or use `woven_damage.ovito`                       |
| Colour coding has no visible effect                       | Applied in a pipeline whose particles and surface are hidden   | Colour the surface vertices (as in `woven_stress.ovito`), or show the particles |
| No `woven_stress.ovito` written                           | Dump made with an old input without the stress columns         | Rerun with the current `in.tensile_sheet`                                   |
| `ModuleNotFoundError: matplotlib`                         | Python without plotting libraries                              | `pip install numpy scipy matplotlib`                                        |

---

## Extending the Model

| Extension                          | What to change                                                                                  |
| ---------------------------------- | ----------------------------------------------------------------------------------------------- |
| Larger or thicker sheet            | `n_rep_x` (length), `n_warp_y` (width), `n_layers` (thickness, 2 or more) in `build_woven.py`  |
| Material properties of the RVE     | Set `periodic = True` and run `in.tensile_periodic` (no grip or edge effects)                   |
| Loading along the weft             | Pull along y: change the grip regions (sheet) or `fix deform` to y and `-pyy` (periodic RVE)   |
| Finer yarns                        | Increase the yarn semi-axes and the spacings `px`, `py` together; keep the bead spacing 1      |
| Different crimp                    | Change `beta` (shape) and `gap`; the script reports any overlap                                |
| Stronger or weaker constituents    | Change `k` and `eps_c` of the bond types in the input                                          |
| Real units                         | Calibrate: set `SIGMA_TO_MPA` so the simulated modulus matches the measured modulus            |
| Beam bonds (bending and twisting)  | `atom_style bpm/sphere`, `bond_style bpm/rotational`, `fix nve/bpm/sphere` (not implemented here) |
| Strain-rate check                  | Halve `erate` and compare the stress-strain curves                                             |

---

## Citation

If you use this code in your research, please cite:

```
@software{mishra_2026_woven_lammps,
  author    = {Mishra, Akshansh},
  title     = {Tensile Testing of a 3D Angle-Interlock Woven Composite Sheet
               Using LAMMPS (BPM Bonded-Particle Model)},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24155083},
  url       = {https://doi.org/10.5281/zenodo.24155083}
}
```

Plain text citation:

> Mishra, A. (2026). *Tensile Testing of a 3D Angle-Interlock Woven Composite Sheet Using LAMMPS (BPM Bonded-Particle Model)* [Computer software]. Zenodo. <https://doi.org/10.5281/zenodo.24155083>

---

## Author

**Akshansh Mishra**
GitHub: <https://github.com/akshansh11>

---

## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />
This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

You are free to:

- **Share**: copy and redistribute the material in any medium or format
- **Adapt**: remix, transform, and build upon the material

Under the following terms:

- **Attribution**: You must give appropriate credit to Akshansh Mishra and provide a link to this repository
- **NonCommercial**: You may not use the material for commercial purposes

Copyright 2026 Akshansh Mishra. All rights reserved for commercial use.
