# NSE Conjugate Heat Transfer Example

This repository contains a deal.II-based 2D flow and temperature solver.
It was originally based on deal.II step-57.
The current setup models incompressible flow in a channel with solid walls and a separate temperature solve over the full fluid-solid domain.

## Build

Configure and build with CMake:

```bash
cmake -S . -B build
cmake --build build
```

The project expects a deal.II installation with `UMFPACK` enabled, as checked in
[`CMakeLists.txt`](CMakeLists.txt).

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
- `--global-refinement` uses uniform global refinement instead of adaptive
  refinement.

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
`heat_exchanger_ref2.vtu`.

## Case Files

Simulation setup is now split into:

- mesh files in [`meshes/`](meshes)
- case files in [`cases/`](cases)

A case file defines:

- which mesh to load
- one subsection per material ID with its kind and properties
- velocity boundary conditions
- temperature Dirichlet boundary conditions

Current examples:

- [`cases/heat_exchanger.prm`](cases/heat_exchanger.prm)
- [`cases/pipe.prm`](cases/pipe.prm)
- [`cases/mms.prm`](cases/mms.prm)
- [`cases/mms_deg2.prm`](cases/mms_deg2.prm)

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

## Manufactured Solution Verification

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
and are wrapped for deal.II in [`src/mms.h`](src/mms.h). MMS mode replaces the
flow and temperature boundary values with the manufactured solution and adds the
corresponding momentum and temperature source terms.

Example MMS runs:

```bash
cd build
./cht_solver --case ../cases/mms.prm --global-refinement --output-dir ../solns/mms_tables
./cht_solver --case ../cases/mms_deg2.prm --global-refinement --output-dir ../solns/mms_tables_deg2
```

When `Use MMS = true`, VTK output includes `velocity_error`, `pressure_error`,
and `temperature_error` fields. The solver also computes convergence data for
`L2_velocity`, `L2_pressure`, and `H1_velocity` after each refinement cycle and
writes both Org and LaTeX tables:

```text
error-<adaptive|global>-q<degree>.org
error-<adaptive|global>-q<degree>.tex
```

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
