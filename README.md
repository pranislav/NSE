# NSE Conjugate Heat Transfer

This repository contains a FEM-based 2D time-independent conjugate heat transfer solver with incompressible flow using deal-II library.
It solves for a velocity field on fluid areas and temperature field on whole mesh.
It was originally based on deal.II step-57.

**Governing equations:**

Incompressible steady-state Navier-Stokes equations:

`(u · ∇)u − ν∆u + ∇p − f = 0`,

`∇ · u = 0`,

where u is the velocity, ν is viscosity, p is pressure and f external body force.

Steady advection-diffusion:

`u · ∇T − κ∆T − s = 0`,

where T is temperature, κ is thermal diffusivity and s is a source term.



## Build

Configure and build with CMake:

```bash
cmake -S . -B build
cmake --build build
```

The project expects a deal.II installation with `UMFPACK` enabled, as checked in
[`CMakeLists.txt`](CMakeLists.txt).

## Build a Portable Bundle

On a Linux machine where deal.II and the other build dependencies are installed,
configure the project and build the `bundle` target:

```bash
cmake -S . -B build
cmake --build build --target bundle
```

This creates:

```text
build/bundle/
  cht_solver
  lib/
```

The `bundle` target copies the solver executable to `build/bundle/`, discovers
its shared-library dependencies with `ldd`, and copies non-system libraries into
`build/bundle/lib/`. This includes deal.II and other non-system runtime
dependencies such as SuiteSparse, BLAS/LAPACK, Fortran runtime libraries, and
C++ runtime libraries when they are resolved as non-system dependencies on the
build machine.

The executable is built with an RPATH containing `$ORIGIN/lib`, so after the
bundle directory is copied to another Linux machine, `cht_solver` searches for
its bundled libraries next to itself in `lib/`.

The bundling step intentionally does not copy standard Linux system libraries
such as `libc`, `libm`, `libpthread`, `librt`, `libdl`, `ld-linux`, and
`linux-vdso`. The target machine therefore still needs a compatible Linux
system and C library, but it does not need deal.II installed system-wide.

Assumptions:

- Bundling is supported on Linux only.
- `ldd` must be available on the build machine.
- `patchelf` is not required; RPATH is set from CMake.
- Input case files and meshes are not embedded in the executable. Copy the
  required files, for example `cases/` and `meshes/`, together with the bundle
  or provide paths to equivalent files on the target machine.

## Run

Run the executable from the `build/` directory and pass a case file:

```bash
cd build
./cht_solver --case ../cases/heat_exchanger.prm
```

Use `--help` to print the available command-line options.

Useful options:

- `--case <file>` selects a case file. Defaults to `../cases/heat_exchanger.prm`.
- `--output-dir <path>` selects where result files are written.
- `--save-mesh` writes the mesh after each refinement cycle.
- `-p`, `--output-partial-solutions` writes intermediate Newton-step output.
- `--global-refinement` uses uniform global refinement instead of adaptive
  refinement.

## Run a Bundle on a Machine Without deal.II

Copy the generated `build/bundle/` directory to the target Linux machine. Also
copy the input files needed by the selected case, usually the corresponding
case file from `cases/` and mesh file from `meshes/`.

Example layout on the target machine:

```text
run/
  bundle/
    cht_solver
    lib/
  cases/
    heat_exchanger.prm
  meshes/
    heat_exchanger.msh
```

Run the bundled executable from `bundle/` or from any working directory:

```bash
cd run/bundle
./cht_solver --case ../cases/heat_exchanger.prm --output-dir ../results/heat_exchanger
```

No `LD_LIBRARY_PATH` setting and no system-wide deal.II installation are needed
for this bundled executable.

If a copied case file contains relative mesh paths, keep the same relative
layout used by the case file or edit the mesh path in the case file before
running.

## Output Directory

You can choose where result files are written with `--output-dir <path>`.
If the directory does not exist, it is created automatically.

Example:

```bash
cd build
./cht_solver --case ../cases/heat_exchanger.prm --output-dir ../results/run_01
```

If `--output-dir` is omitted, output files are written to the current working directory.
Mesh outputs saved with `--save-mesh` are written after each refinement and are
named from the input mesh stem plus the refinement count, for example
`heat_exchanger_ref2.vtk`.

## Case Files

Simulation setup is now split into:

- mesh files in [`meshes/`](meshes)
- case files in [`cases/`](cases)

A case file defines:

- solver parameters (Reynolds number, polynomial degree of basis functions, number of refinements and more)
- which mesh to load
- material type (fluid/solid) and thermal diffusivity of defined mesh areas
- velocity boundary conditions
- temperature boundary conditions

Current examples:

- [`cases/pipe.prm`](cases/pipe.prm) - simple pipe with solid walls
- [`cases/heat_exchanger.prm`](cases/heat_exchanger.prm) - two parallel tubes with flow in opposite direction, thermally connected, with one hot inflow and one cold inflow
- [`experiment_layout.prm`](cases/experiment_layout.prm) - a particular experiment geometry
- [`cases/mms.prm`](cases/mms.prm) - case for code verification by method of manufactured solutions (MMS)
- [`cases/mms_deg<n>.prm`](cases/mms_deg2.prm) where *n* in (2, 3) - MMS of higher degree

## Case File Format

Case files use deal.II parameter syntax and must define:

```prm
subsection Mesh
  set File = ../meshes/<mesh_name>.msh
end

subsection Solver
  set Reynolds number = <reynolds_number>
  set Gamma = <gamma>
  set Polynomial degree = <degree>
  set Refinement cycles = <refinement_cycles>
  set Use MMS = <true|false>
end

subsection Materials
  set Ids = <comma_separated_material_ids>

  subsection Material <id>
    set Kind = <fluid|solid>
    set Thermal diffusivity = <thermal_diffusivity>
  end
end

subsection Boundary conditions
  set Velocity Dirichlet = <boundary_id:type:component:coordinate:value; ...>
  set Temperature Dirichlet = <boundary_id:value; ...>
end
```

Velocity boundary types:

- `constant`
- `parabolic`

Velocity entry format:

- `boundary_id:type:component:coordinate:value`
- `component`: `0` for `u_x`, `1` for `u_y`
- `coordinate`: axis used by the profile, `0` for `x`, `1` for `y`
- `value`: imposed constant value or parabola peak value

No-slip boundaries are represented as zero-valued constant velocity Dirichlet
entries, one per constrained velocity component.

## Code Verification by Method of Manufactured Solutions

The solver has a method of manufactured solutions (MMS) mode for verifying the
Navier-Stokes solver on the unit square (extension of MMS for the temperature part is in progress TODO). Enable it from a
case file with:

```prm
subsection Solver
  set Use MMS = true
end
```

The current manufactured fields are defined in [`scripts/expressions_mms.py`](scripts/expressions_mms.py):

- `u = (sin(pi x) cos(pi y), -cos(pi x) sin(pi y))`
- `p = cos(pi x) sin(pi y)`
- `T = cos(2 pi x) sin(pi y)`

The generated C++ functions live in [`src/mms_generated.h`](src/mms_generated.h)
and are wrapped for deal.II in [`src/mms.h`](src/mms.h).
Note that in case of changing the manufactured solution expressions it is necessary to also regenerate the `mms_generated.h` by running `python3 scripts/generate_mms.py`.
MMS mode defines the flow and temperature boundary values in accordance with the manufactured solution (so any user-defined boundary conditions in case file will be overrode TODO verify experimentally). It also adds the corresponding momentum and temperature source terms in the assembly.
It is recommended to run MMS mode with global refinement option so the convergence rates are more straightforward to compute.

Example MMS run:

```bash
cd build
./cht_solver --case ../cases/mms.prm --global-refinement --output-dir ../solns/mms
```

When `Use MMS = true`, VTK output includes `velocity_error`,
`pressure_error`, and `temperature_error` fields. The solver also computes
convergence data for `L2_velocity`, `L2_pressure`, and `H1_velocity` after each
refinement cycle and writes both Org and LaTeX tables.


```text
error-<adaptive|global>-q<degree>.org
error-<adaptive|global>-q<degree>.tex
```

By default, the solver writes only the final converged solution from the last
refinement cycle. Use `-p` or `--output-partial-solutions` to write the
intermediate Newton-step outputs.

Generate or refresh the MMS expressions with:

```bash
python3 scripts/generate_mms.py
```

This updates [`src/mms_generated.h`](src/mms_generated.h) and
[`mms_report.tex`](mms_report.tex). Plot convergence rates from an Org table
with:

```bash
python3 scripts/convergence_rate_analysis.py solns/mms_tables/error-global-q1.org
```

The plotting script writes one PNG per error column next to the input table.
Note that the script assumes global refinement option was used.
