# LAMMPS Workflow: Compression → Shear-Rate Sweep → Rest (Gay–Berne Ellipsoids)

This repository contains a 3-stage LAMMPS workflow for a **polydisperse mixture of Gay–Berne (GB) ellipsoids** (aspherical particles) meant to mimic qualitative aspects of layered/anisotropic materials (e.g., smectite-like ordering). The workflow:

1) **Build & compress** an initially dilute system to a target pressure.  
2) **Shear** the compressed state at controlled shear rates (sweep) while monitoring stress, alignment, and plasticity.  
3) **Rest/relax** the sheared structure and observe structural recovery.

In addition to the original mechanics (NVT → NPH compression → NPH+deform shear → relaxation), this repo adds **novel analysis** that is not typically included in “standard” input decks:

- **Nematic alignment order parameter** \(S(t)=⟨P2(\cos\theta)⟩\) and **spatial profile** \(S(z)\)
- **Non-affine displacement** \(⟨D^2_{NA}⟩\) under shear (plasticity proxy)
- **Automated shear-rate sweep** with consistent naming of outputs

---

## Files

- `in.1_compression`  
  Builds the polydisperse GB system, equilibrates, then compresses through a pressure schedule and writes restart files.

- `in.2_shear_sweep_press1000`  
  Reads `mix_press1000.restart`, converts the box to triclinic, then shears at several strain rates \(\dot\gamma\) while controlling normal stress (\(p_{zz}\)). Outputs stress, strain, alignment, non-affine motion, and velocity profiles.

- `in.3_rest_press1000`  
  Reads a chosen sheared restart (default: the highest-rate tag) and relaxes at constant volume (NVE + thermostat). Optionally relax at constant normal stress with NPH.

---

## Requirements

LAMMPS built with:

- **ASPHERE** package (for ellipsoid integration: `atom_style ellipsoid`, `nve/asphere`, `nph/asphere`)
- **OPENMP** package (for `gayberne/omp`), optional but recommended

The original README for these inputs (provided with the files you started from) recommends **LAMMPS 23Jun2022** and an OpenMP build.

---

## Physical model (what is being simulated)

### Particle model: Gay–Berne ellipsoids
Each particle is an **aspherical ellipsoid** interacting via the **Gay–Berne potential**, which is an orientation-dependent generalization of Lennard–Jones. It favors certain alignments and can produce **liquid crystalline ordering** (nematic/smectic-like) depending on aspect ratio, density, and temperature.

### Polydispersity
The initial configuration contains **three particle types** with different masses and shapes:

- Type 1: shape \((a,b,c)=(1.0,1.0,0.2)\), mass 1.0
- Type 2: shape \((1.2,1.2,0.2)\), mass 1.44
- Type 3: shape \((1.4,1.4,0.2)\), mass 1.96

This polydispersity is introduced to prevent crystallization and to encourage glassy/disordered packing at high pressure.

### Units
All scripts use `units lj` (reduced LJ units). That means:

- Temperature, energy, time, pressure are all in **reduced units**.
- A timestep of `0.00025` is a reduced time step.

---

## Stage 1 — Compression (`in.1_compression`)

### Goal
Create a well-defined packed structure by equilibrating at constant volume, then compressing through a sequence of pressures, writing a restart at each pressure.

### Ensemble sequence
1. **NVT-like equilibration** implemented as:
   - `fix nve/asphere` (integrator)
   - `fix langevin ... angmom 1.0` (thermostat including rotational DOF)

2. **Compression using NPH**:
   - `fix nph/asphere` ramps pressure from one target to the next
   - then equilibrates at constant target pressure

### Temperature and pressure definitions
Aspherical particles have both translational and rotational kinetic energy. The scripts define:

- `compute temp_trans` : translational temperature
- `compute temp/asphere` : combined rotational+translational temperature
- Pressure uses `pressure temp_trans` (common choice: define pressure from translational motion to avoid mixing thermostatting of rotations into the pressure compute).

### Output: Alignment order
This stage now outputs:

- **Time series**: `compression.S_time.dat`
- **Profile**: `compression_novel.S_profile.dat`

These quantify how compression induces orientational ordering.

---

## Stage 2 — Shear-rate sweep (`in.2_shear_sweep_press1000`)

### Goal
Starting from a compressed restart at `press = 1000`, impose steady simple shear at different strain rates \(\dot\gamma\). Measure stress response, velocity profile, ordering, and non-affine motion.

### Triclinic box and shear deformation
To shear a periodic cell in LAMMPS, the simulation box is converted to **triclinic**, then the **tilt factor** is driven:

- `change_box all triclinic ...`
- `fix deform 1 xz erate ${gdot} remap v`

This imposes a time-dependent tilt `xz(t)` corresponding to simple shear in the x–z plane.

The instantaneous (engineering) shear strain is tracked as:

- \(\gamma(t) = xz/lz\)

### Ensembles under shear
- Pressure is controlled only in z:
  - `fix nph/asphere z 1000 1000 0.25`

- Temperature control uses `temp/deform`:
  - under shear there is a streaming velocity profile; `temp/deform` subtracts it so the thermostat does not artificially remove flow energy.

### Shear-rate sweep
Instead of a single \(\dot\gamma\), the script runs a list:

- 0.001, 0.003, 0.01, 0.03, 0.04

Each run writes separate outputs tagged by `${tag}` (001, 003, 010, 030, 040).

### Output A: Nematic alignment \(S(t)\) and \(S(z)\)
We compute a unit vector (body axis) from the quaternion and form:

- \(P_2 = \frac{1}{2}(3u_z^2 - 1)\)
- \(S(t)=\langle P_2\rangle\)

Outputs:

- `mix4000_shear${tag}.S_time.dat`  (time series)
- `mix4000_shear${tag}.S_profile.dat` (profile vs z)

Interpretation:

- \(S \approx 0\) : isotropic orientations
- \(S > 0\) : alignment with z-axis (nematic-like)
- \(S < 0\) : alignment perpendicular to z-axis

### Novel output B: Non-affine displacement \(⟨D^2_{NA}⟩\)
Under affine shear, positions would follow:

- \(x_{\mathrm{aff}}(t) = x_0 + \gamma(t) z_0\)

Non-affine displacement removes this part:

- \(\Delta x_{NA} = x - x_0 - \gamma z_0\)
- \(D^2_{NA} = \Delta x_{NA}^2 + \Delta y^2 + \Delta z^2\)

We output:

- `mix4000_shear${tag}.nonaffine.dat` containing \(\gamma\) and \(⟨D^2_{NA}⟩\)

Interpretation:

- Small \(⟨D^2_{NA}⟩\): mostly elastic/affine response
- Rapid growth: plastic rearrangements / yielding / shear band activity

### Standard outputs (also kept)
- Stress tensor components (especially `pxz` shear stress)
- Normal stress `pzz`
- Velocity profile `mix4000_shear${tag}.profile`
- Dumps with positions and quaternions

---

## Stage 3 — Rest / Relaxation (`in.3_rest_press1000_novel`)

### Goal
After shear, stop deforming and let the system relax. This probes structural recovery, stress relaxation, and possible re-ordering.

### Default rest condition
- Constant volume with `fix nve/asphere` + Langevin thermostat

### Optional rest at constant normal stress
Uncomment the NPH fix:

```lammps
#fix 4 all nph/asphere z 1000.0 1000.0 0.25
#fix_modify 4 press pres_trans
