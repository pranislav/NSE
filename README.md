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

The project expects a deal.II installation with `UMFPACK` enabled, as checked in [`CMakeLists.txt`](/home/branislav/Vysoka/phd?/deal/NSE/CMakeLists.txt).

## Run

Run the executable from the `build/` directory and pass a case file:

```bash
cd build
./cht_solver --case ../cases/heat_exchanger.prm
```

Use `--help` to print the available command-line options.

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

- mesh files in [`meshes/`](/home/branislav/Vysoka/phd?/deal/NSE/meshes)
- case files in [`cases/`](/home/branislav/Vysoka/phd?/deal/NSE/cases)

A case file defines:

- which mesh to load
- one subsection per material ID with its kind and properties
- velocity boundary conditions
- temperature Dirichlet boundary conditions

Current examples:

- [`cases/heat_exchanger.prm`](/home/branislav/Vysoka/phd?/deal/NSE/cases/heat_exchanger.prm)
- [`cases/pipe.prm`](/home/branislav/Vysoka/phd?/deal/NSE/cases/pipe.prm)

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
  set Adaptive refinement cycles = <adaptive_refinement_cycles>
end

subsection Materials
  set Ids = <comma_separated_material_ids>

  subsection Material <id>
    set Kind = <fluid|solid>
    set Thermal conductivity = <thermal_conductivity>
  end
end

subsection Boundary conditions
  set No-slip boundaries = <comma_separated_boundary_ids>
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
