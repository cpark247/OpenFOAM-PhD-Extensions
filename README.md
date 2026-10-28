# OpenFOAM-PhD-Extensions

Custom solvers, turbulence and population-balance models, boundary conditions and utilities that I developed during my doctoral research on gas–liquid flows in hydraulic machinery and diverging channels at the Chair of Hydraulic Fluid Machinery, Ruhr University Bochum. Everything targets **OpenFOAM Foundation v9**.

| | |
|---|---|
| **Used in** | C. Park, M. Mansour, D. Thévenin, R. Skoda: *Passive Mitigation of Gas Accumulation by Surface Structures in Diverging Flow Channels: A Three-Dimensional Two-Phase Flow Simulation Study*, [SSRN preprint](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7649860) |
| **Validation cases** | [turbulence-validation-OF9](https://github.com/cpark247/turbulence-validation-OF9): channel, square duct and asymmetric-diffuser benchmarks run with these tools |
| **Paper companion** | [gas-accumulation-diverging-channel](https://github.com/cpark247/gas-accumulation-diverging-channel) |
| **Licence** | GPL-3.0 (derived from OpenFOAM) |

---

## Contents

### Multiphase: Euler–Euler with population balance (`bubbleqmom-dev/`)

A hybrid multiphase solver family built on `multiphaseEulerFoam`, extended for **scale-resolving** (LES, SST-SAS) as well as RANS simulation of bubbly flows with population balance modelling.

| Component | What it does |
|---|---|
| `multiphaseEulerFoamHSM_RANS` | RANS variant with run-time time-averaging of the flow fields, per-size-class **coalescence and breakup rate** storage and averaging, and Sauter mean diameter (d₃₂) output |
| `multiphaseEulerFoamHSM_calcResLES` | Same for LES: resolved, modelled and total turbulence statistics plus per-size-class population-balance rates and d₃₂ |
| `multiphaseEulerFoamHSM_calcResLES_mono` | Resolved-field statistics without the population-balance rate fields; supports LES (k from ν_t and Δ) and SST-SAS (k and ε from the RANS fields) |
| `LuoSvendsenForSAS2`, `PrinceBlanchForSAS2` | Luo–Svendsen breakup and Prince–Blanch coalescence kernels that use the **total** (resolved + modelled) dissipation rate, so the kernels remain consistent in scale-resolving runs where the modelled part alone is too small |
| `LuoSvendsenForSST`, `PrinceBlanchForSST` | The same kernels driven by the RANS-modelled dissipation rate |
| `phaseModelledTKE` | Function object for the modelled turbulence kinetic energy / dissipation of a phase, with time averaging |

### Turbulence models (`turbModels_mod/`)

Compiled into `libincompressibleMomentumTransportModels_mod`; load with `libs ("libincompressibleMomentumTransportModels_mod.so");` in `controlDict`.

| Model | Change relative to OpenFOAM-9 |
|---|---|
| `kOmegaSST_ext` | Standard k-ω SST that additionally writes its production term `Pk_` for analysis |
| `kOmegaSSTSAS_ext` | k-ω SST-SAS (Egorov & Menter 2008) with a configurable exponent `nExponent` on the length-scale ratio L/L_vK in the Q_SAS source (default 2 = original model), and output of `L_`, `Lvk_`, `L_by_Lvk_` and `Qsas_` to diagnose when and where the model switches to scale-resolving behaviour |
| `kOmegaSSTIDDES`, `BSL_EARSM` | Work in progress: source included, **not registered/compiled** by default |

### Single-phase solvers with run-time statistics

Run-time time-averaging is done inside the solver rather than via the `fieldAverage` function object, to work around limitations I ran into with the function-object route for these quantities.

| Solver | Additional output |
|---|---|
| `pimpleFoam_calcMean` | Time-averaged p and U, modelled k, ε and ν_t, and the resolved Reynolds stress tensor; for URANS and hybrid runs |
| `pimpleFoam_calcRes` | LES-specific: resolved, modelled and total k and ε, individual Reynolds-stress components (uu … vw), and Δ/η (filter width over Kolmogorov scale) as a mesh-resolution indicator |
| `pimpleFoam_calcRes_S2W2` | As `pimpleFoam_calcMean`, plus time-averaged strain and rotation invariants S², W² and S²/W² for SST-SAS analysis |

### Boundary conditions

| Boundary condition | What it does |
|---|---|
| `porousBafflePressureAMI` | `porousBafflePressure` re-implemented on `fixedJumpAMI`, so it works on **non-conformal** `cyclicAMI` / `cyclicACMI` interfaces. Models a hydraulic resistance (Darcy–Forchheimer pressure jump) on a pair of periodic patches without meshing the physical geometry. This is the core of the *Virtual Hole Method* for holes in centrifugal-pump impeller blades: the hole's effect is set in a dictionary instead of meshing each hole geometry, which cut design iteration cycles from weeks to minutes. Works with single-phase solvers and `multiphaseEulerFoam` |
| `periodicTimeVaryingMappedFixedValue` | Like `timeVaryingMappedFixedValue`, but loops back to the first time once the end of the supplied boundary data is reached. Used to recycle a finite precursor inflow database for long LES runs |

### Utilities

| Utility | What it does |
|---|---|
| `perturbUChannel`, `perturbUCylinder`, `perturbUGeneric` | Velocity initialisation with near-wall streaks to trigger transition in LES (de Villiers 2006; see attribution) |
| `perturbUDuct`, `perturbURectDuct` | New variants for square and rectangular ducts with arbitrary duct-centre position and geometry checks |
| `calcEpsilon`, `calcOmega` | Convert between ω and ε (ε = C_μ k ω) when initialising or switching turbulence models |
| `calcRST` | Adds resolved and modelled parts into the total Reynolds stress tensor |

---

## Building

Requires OpenFOAM Foundation v9 on Linux (developed on HPC clusters with the Intel and GCC toolchains).

```bash
source /opt/openfoam9/etc/bashrc     # or your own OpenFOAM-9 installation
git clone https://github.com/cpark247/OpenFOAM-PhD-Extensions.git
cd OpenFOAM-PhD-Extensions
./Allwmake -j
```

Libraries are installed to `$FOAM_USER_LIBBIN` and executables to `$FOAM_USER_APPBIN`. `./Allwclean` removes the build output. Components can also be built individually with `wmake <directory>`.

Example use of the boundary condition in `0/p` (both patches of the cyclicACMI pair use it; add `libs ("libporousBafflePressureAMI.so");` to `system/controlDict`):

```c++
hole_inlet
{
    type            porousBafflePressureAMI;
    patchType       cyclicACMI;
    jump            uniform 0;
    D               0;          // Darcy coefficient    [1/m^2]; 0 is fine when inertia dominates
    I               100;        // inertial coefficient [1/m]
    length          0.015;      // porous length in flow direction [m]
    // multiphaseEulerFoam only:
    // phi  phi.water;  rho  rho.water;  phase  water;
    // optional: positive false; jumpMax 30;
    value           uniform 0;
}
```

The full option table is in the header of `porousBafflePressureAMIFvPatchField.H`.

---

## Attribution

This repository contains modified and redistributed third-party GPL code. My own work is the set of extensions listed above; everything else belongs to its original authors.

- **OpenFOAM** – © OpenFOAM Foundation, GPL-3.0. Most files start from OpenFOAM-9 sources; headers keep the original copyright.
- **Quadrature-based moment methods** (`bubbleqmom-dev/src/quadratureMethods`, `Vandermonde`, `eigenSolver`, `mappedList`, …) – from [OpenQBMM](https://www.openqbmm.org/) by A. Passalacqua and co-authors, GPL-3.0.
- **Hybrid multiphase solver base** (`multiphaseEulerFoamHSM`) – Hundshagen (2023), Chair of Hydraulic Fluid Machinery, Ruhr University Bochum. My variants are listed above.
- **perturbU** – E. de Villiers, *The Potential of Large Eddy Simulation for the Modeling of Wall Bounded Flows*, PhD thesis, Imperial College London, 2006; maintained at [wyldckat/perturbU](https://github.com/wyldckat/perturbU) (GPL-2.0, see `perturbU/LICENSE`). The duct variants are mine.
- Model references are given in the header of each source file.

## Citation

If you use this code, please cite the paper above and this repository (see `CITATION.cff`).

## Licence

GPL-3.0, see [`LICENSE`](LICENSE). This offering is not approved or endorsed by the OpenFOAM Foundation or OpenCFD Ltd., producer and distributor of the OpenFOAM software and owner of the OpenFOAM® trademark.
