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

## Output Directory

You can choose where result files are written with `--output-dir <path>`.
If the directory does not exist, it is created automatically.

Example:

```bash
cd build
./cht_solver --case ../cases/heat_exchanger.prm --output-dir ../results/run_01
```

If `--output-dir` is omitted, output files are written to the current working directory.

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
