[![DOI](https://zenodo.org/badge/1281491550.svg)](https://doi.org/10.5281/zenodo.20939259)

# MegathrustRS

MATLAB code for one-dimensional rate-and-state friction simulations on a discretized subduction-interface fault. This code is based largely on original matlab codes that became the framework for FDRA (Fault Dynamics with a Radiation Damping Approximation): Segall, P., & Bradley, A. M. (2012). The role of thermal pressurization and dilatancy in controlling the rate of fault slip. Journal of Applied Mechanics, 79(3), 031013. https://doi.org/10.1115/1.4005896.

The code builds an elastic stress-interaction matrix for a prescribed fault geometry, compresses the interaction kernel with H-matrix matrix-vector products, defines depth-dependent frictional and hydraulic properties, runs quasi-dynamic slip simulations, and provides plotting scripts for the resulting slip, stress, state, and pore-pressure evolution.

This repository is organized around a four-step workflow:

1. Build the fault geometry and elastic Green's functions.
2. Define the friction, pore-pressure, dilatancy, and loading parameters.
3. Run the rate-and-state simulation.
4. Load and plot the saved outputs.

An optional fifth step allows an external shear-stress perturbation to be imposed on a previously computed simulation state.

## Contents

```text
RS_friction_subduction/
├── Brief_Instructions.txt                  # Short workflow notes
├── setup_geometry_GFs.m                    # Build fault geometry and elastic stress Green's functions
├── setup_simulation.m                      # Main simulation setup file
├── run_simulation.m                        # Run rate-and-state/dilatancy ODE simulation
├── load_and_plot_outputs.m                 # Read and plot saved output files
├── run_external_stress_change.m            # Restart from a previous run with imposed stress change
├── run_simulation_large_Linf.m             # Alternative simulation driver for large-Linf setup
├── setup_simulation_large_Linf.m           # Alternative simulation setup
├── sample_setup_scripts/                   # Example setup files for common model configurations
│   ├── setup_simulation_simple_sequence.m
│   ├── setup_simulation_ikari.m
│   ├── setup_simulation_dilatancy1.m
│   └── setup_simulation_dilatancy2.m
├── tools/                                  # MATLAB helper functions for geometry, stress, ODEs, plotting
├── hmmvp0.16/                              # HMMVP package for Mac/Linux-style builds
├── hmmvp0.16_win/                          # HMMVP package configured for Windows builds
└── lib/                                    # Windows libraries used by HMMVP compilation
```

## Requirements

- MATLAB.
- A working MATLAB MEX compiler for C++ code if the HMMVP MEX files need to be compiled on your system.
- The bundled HMMVP package, or an already compiled compatible HMMVP installation.

The code uses the HMMVP package by Andrew M. Bradley to accelerate matrix-vector products with a hierarchical-matrix approximation to the elastic stress kernel. HMMVP is distributed under the Eclipse Public License 1.0; see `hmmvp0.16/readme.txt` for installation notes and citation information.

## Quick start

From MATLAB, start in the repository root directory.

```matlab
% 1. Build the geometry and elastic Green's functions.
setup_geometry_GFs

% 2. Define the frictional/hydraulic model and initial conditions.
setup_simulation

% 3. Run the simulation.
run_simulation

% 4. Plot outputs.
load_and_plot_outputs
```

The geometry/Green's-function setup only needs to be rerun when the fault geometry, patch length, or H-matrix tolerance changes. For parameter studies using the same geometry, typically rerun only `setup_simulation.m` and `run_simulation.m`.

## Installation and setup

### 1. Clone or download the repository

```bash
git clone <repository-url>
cd RS_friction_subduction
```

### 2. Compile HMMVP if needed

The repository includes HMMVP source code and some precompiled MEX files. These compiled files may not work on every operating system, MATLAB version, or processor architecture. If MATLAB cannot find or execute `hm_mvp` or `hmcc_hd_mex`, compile HMMVP locally.

For Mac/Linux-style systems:

```matlab
cd hmmvp0.16
make
cd ..
```

For Windows-style builds:

```matlab
cd hmmvp0.16_win
make
cd ..
```

If compilation fails, check that MATLAB has a configured C++ compiler:

```matlab
mex -setup C++
```

Then retry the HMMVP `make` script. Additional HMMVP installation notes are in `hmmvp0.16/readme.txt`.

### 3. Confirm MATLAB paths

Most top-level scripts add the required paths internally. In particular:

- `setup_geometry_GFs.m` adds `tools/` and `hmmvp0.16/`.
- `setup_simulation.m` adds `tools/` and selects `hmmvp0.16_win/` on Windows or `hmmvp0.16/` otherwise.

If you move directories or use an external HMMVP installation, update the relevant `addpath` statements.

## Workflow in detail

### Step 1: Build geometry and elastic Green's functions

Run:

```matlab
setup_geometry_GFs
```

This script defines the subduction-interface geometry, discretizes it into fault patches, and builds an elastic shear-stress interaction kernel. The interaction kernel is stored as an H-matrix so that later simulations can use fast matrix-vector products rather than explicitly storing or multiplying a dense matrix.

Key input parameters near the top of `setup_geometry_GFs.m` are:

| Parameter | Meaning | Units |
|---|---|---|
| `faultdip_trench` | Dip of the interface at the trench | degrees |
| `faultdip_bottom` | Dip of the interface at the bottom of the modeled interface | degrees |
| `x_trench` | Horizontal position of the trench | km |
| `x_bottom` | Horizontal position of the bottom of the modeled interface | km |
| `z_bot` | Depth of the bottom of the modeled interface | km |
| `pL` | Maximum patch length used for elastic calculations | km |
| `name_elastic` | Output `.mat` filename for geometry and H-matrix metadata | string |

The default output file is:

```text
elastic_GFs.mat
```

The saved file contains, among other variables:

- `FaultPatches` — patch endpoint coordinates, `[top_x, top_z, bottom_x, bottom_z]`.
- `Geometry` — structure containing the geometry settings.
- `gamb` — structure containing H-matrix and ODE-output metadata.

The script also writes H-matrix data files such as:

```text
shear_Hmat_ne<hmat-settings>.dat
```

These files are needed by later simulations. Do not delete them unless you plan to rebuild the Green's functions.

### Step 2: Define the frictional and hydraulic model

Run:

```matlab
setup_simulation
```

This script loads the elastic Green's-function file, defines the long-term loading rate and frictional/hydraulic parameters, interpolates those properties onto the fault patches, computes derived quantities such as effective normal stress and nucleation length scale, and generates diagnostic plots.

Most user edits should be made in the `INPUT section` at the top of `setup_simulation.m`.

Important inputs include:

| Parameter | Meaning | Units |
|---|---|---|
| `elastic_name` | Name of the `.mat` file created by `setup_geometry_GFs.m` | string |
| `savename` | Output filename for the simulation time series | string |
| `srate` | Long-term slip/loading rate | m/yr |
| `mu_0` | Nominal coefficient of friction | dimensionless |
| `param_depth` | Depth control points for parameter interpolation | km |
| `a`, `b` | Rate-and-state friction parameters at `param_depth` | dimensionless |
| `d_c` | Critical slip distance at `param_depth` | m |
| `pore_pressure_fraction` | Pore pressure as a fraction of lithostatic pressure | dimensionless |
| `max_eff` | Maximum allowed effective normal stress | Pa |
| `c_star` | Pore-pressure diffusivity parameter | 1/yr |
| `beta` | Bulk compressibility | 1/Pa |
| `eps` | Dilatancy coefficient | dimensionless |
| `stopt` | Simulation stop time | yr |

Depth in this setup is depth below the seafloor/ground surface, not depth below sea level. The first and last `param_depth` values should correspond to the top and bottom edges of the modeled interface. The default bottom value is set from the geometry:

```matlab
param_depth(end) = -FaultPatches(end,4);
```

#### Fixed versus velocity-dependent friction

The code supports two ways to define `a-b`:

1. Fixed `a` and `b` values assigned as functions of depth.
2. A velocity-dependent `a-b` formulation based on Ikari-style behavior for selected patches.

The logical vector `use_ikari` determines which depth intervals use the velocity-dependent formulation:

```matlab
use_ikari(i) = true;   % use velocity-dependent a-b in interval i
use_ikari(i) = false;  % use fixed a and b in interval i
```

When `use_ikari=true`, the relevant parameters are:

| Parameter | Meaning |
|---|---|
| `a_b_slope` | Slope of `a-b` versus `log10` slip rate |
| `a_b_cut` | Minimum allowed `a-b` |
| `logVcut` | Slip rate cutoff used in the velocity-dependent relation |

#### Diagnostic plots

`setup_simulation.m` creates plots showing:

- `a`, `b`, and velocity-dependent bounds versus depth.
- `a-b` on the fault interface.
- Lithostatic stress, pore pressure, and effective normal stress versus depth.
- Nucleation length scale `h*` compared with model patch length.
- Velocity-dependent `a`, `b`, and `a-b` behavior where `use_ikari=true`.

Check the `h*` plot before running long simulations. Patch lengths should be smaller than the relevant nucleation length scale.

### Step 3: Run the simulation

Run:

```matlab
run_simulation
```

This script initializes the H-matrix matrix-vector product, computes the elastic stressing rate associated with steady loading, and integrates the rate-and-state equations with dilatancy using:

```matlab
get_ratestate_dilatancy(...)
```

The simulation can take minutes to hours depending on the number of patches, frictional properties, and integration behavior. The ODE output function saves results incrementally to disk, so a stopped or interrupted simulation will usually retain output up to the last saved step.

The output filename is controlled by `savename` in `setup_simulation.m` and by:

```matlab
gamb.ode.savefn = [dr '/' savename];
```

in `run_simulation.m`.

The parameter:

```matlab
gamb.ode.saveEvery = 10;
```

controls how frequently solver steps are saved.

### Step 4: Load and plot outputs

Before plotting, make sure the simulation setup variables are available in the MATLAB workspace. If needed, rerun:

```matlab
setup_simulation
```

Then edit the output filename at the top of `load_and_plot_outputs.m`:

```matlab
savefn = 'ode_elastic_simple_sequence';
```

and run:

```matlab
load_and_plot_outputs
```

The script reads the saved time series with:

```matlab
x = SaveStreamData('Read',savefn,stride,[],1);
```

where `stride` controls how many saved time steps are skipped while loading. Increasing `stride` can make plotting faster for large output files.

The saved state vector is unpacked into:

| Variable | Meaning |
|---|---|
| `vel` | Slip rate on each patch |
| `theta` | State variable |
| `tau` | Shear stress |
| `phi` | Auxiliary dilatancy/state variable |
| `p` | Pore pressure |
| `B` | Depth- and/or velocity-dependent `b` parameter |
| `slip` | Cumulative slip |

The plotting script produces time-depth or time-distance plots of:

- `b` value.
- `log10` slip rate.
- Cumulative slip.
- Shear stress.
- State variable.
- Pore pressure.
- Cumulative slip evolution.
- Slip-rate time series at selected distances along the fault.

## Optional workflow: impose an external stress change

The script `run_external_stress_change.m` restarts a simulation from an existing output file and applies a prescribed shear-stress perturbation.

Before running it, first run the setup file corresponding to the original simulation so that geometry and model variables are available in the workspace. Then edit the input section of `run_external_stress_change.m`.

Important inputs are:

| Parameter | Meaning | Units |
|---|---|---|
| `savefn_existing` | Existing ODE output file used to extract initial conditions | string |
| `savename` | Output filename for the restarted simulation | string |
| `t_impose` | Time in the original simulation at which the stress perturbation is imposed | yr |
| `delta_tau` | Shear-stress change on each fault patch | Pa |
| `duration` | Ramp duration for applying the stress change | yr |
| `stopt` | Duration of restarted simulation | yr |

Example stress perturbation:

```matlab
delta_tau = 10^5*exp(-(depths-20000).^2/5000^2);
```

This imposes a Gaussian-shaped shear-stress change centered near 20 km depth.

## Example setup scripts

The directory `sample_setup_scripts/` contains setup files that can be copied to the repository root or used as templates for different model configurations.

| File | Purpose |
|---|---|
| `setup_simulation_simple_sequence.m` | Simple depth-dependent rate-and-state setup for long sequence simulations |
| `setup_simulation_ikari.m` | Example using velocity-dependent Ikari-style `a-b` behavior in selected depth intervals |
| `setup_simulation_dilatancy1.m` | Dilatancy example with a shallow/intermediate dilatant interval |
| `setup_simulation_dilatancy2.m` | Dilatancy example with a deeper dilatant interval and longer interface |

A typical way to use one of these examples is:

```matlab
copyfile('sample_setup_scripts/setup_simulation_ikari.m','setup_simulation.m')
setup_simulation
run_simulation
```

Alternatively, run the sample setup script directly, then run `run_simulation.m`, provided the expected workspace variables have been created.

## Units and sign conventions

- Distances in geometry setup are in kilometers.
- Patch lengths used internally for frictional calculations are converted to meters.
- Slip rates are in meters per year.
- Simulation time is in years.
- Stresses and pore pressures are in pascals.
- Depth is treated as positive downward for most diagnostic plotting, although the fault-patch coordinates use negative `z` values for depth in the geometry arrays.
- `sigma` is lithostatic normal stress.
- `pinf` is the nominal pore pressure outside the deforming fault zone.
- Effective normal stress is `sigma - pinf`.

## Main output files

The workflow produces three main classes of output:

1. Geometry and Green's-function metadata:

   ```text
   elastic_GFs.mat
   ```

2. H-matrix data files used for fast matrix-vector products:

   ```text
   shear_Hmat_*.dat
   normal_Hmat_*.dat   # only if normal-stress H-matrix generation is enabled
   ```

3. ODE output files with the saved simulation state:

   ```text
   ode_elastic_*
   ```

The ODE files are read with `SaveStreamData('Read',...)`, not with MATLAB's standard `load` command.

## Editing guide for new simulations

For most new experiments, edit only the input sections of these files:

1. `setup_geometry_GFs.m`
   - Change the interface geometry and patch length.
   - Rerun this only when geometry or discretization changes.

2. `setup_simulation.m`
   - Change frictional properties, pore-pressure structure, dilatancy parameters, loading rate, simulation duration, and output filename.
   - Rerun this whenever model parameters change.

3. `run_simulation.m`
   - Usually only change `gamb.ode.saveEvery`.
   - The rest of the file normally does not need to be edited.

4. `load_and_plot_outputs.m`
   - Change `savefn`, `stride`, and any plot-specific settings.

5. `run_external_stress_change.m`
   - Change the source output file, perturbation time, stress-change vector, ramp duration, and restarted output filename.

## Troubleshooting

### MATLAB cannot find `hm_mvp`

Check that the correct HMMVP directory has been added to the MATLAB path:

```matlab
addpath hmmvp0.16       % Mac/Linux-style directory
% or
addpath hmmvp0.16_win   % Windows-style directory
```

If the function is still unavailable, compile HMMVP using its `make.m` file.

### MEX file has the wrong architecture

Precompiled `.mex*` files are platform-specific. Recompile HMMVP on the machine and MATLAB version where the simulations will be run.

### Simulation is very slow

Possible causes include:

- Very small patch length `pL`, which increases the number of elements.
- Frictional parameters that produce very small nucleation scales or rapid slip transients.
- Very strict H-matrix tolerance.
- Very frequent output saving through `gamb.ode.saveEvery`.

Start with a coarser geometry or shorter `stopt` while testing a new model.

### Output file is too large

Increase `gamb.ode.saveEvery` in `run_simulation.m`, or increase `stride` in `load_and_plot_outputs.m` when reading outputs for plotting.

### Plots fail after restarting MATLAB

`load_and_plot_outputs.m` requires variables created by `setup_simulation.m`, including geometry, patch lengths, depths, and loading parameters. Rerun the setup script before loading outputs.

### Patch length is larger than `h*`

The setup script plots nucleation length scale `h*` and patch length versus depth. If patch length is larger than `h*` in the velocity-weakening region, decrease `pL` in `setup_geometry_GFs.m`, rebuild the Green's functions, rerun `setup_simulation.m`, and rerun the simulation.

