# NSE Conjugate Heat Transfer Example

This repository contains a deal.II-based 2D flow and temperature solver built from `step57base.cc`.
The current setup models incompressible flow in a channel with solid walls and a separate temperature solve over the full fluid-solid domain.

## Build

Configure and build with CMake:

```bash
cmake -S . -B build
cmake --build build
```

The project expects a deal.II installation with `UMFPACK` enabled, as checked in [`CMakeLists.txt`](/home/branislav/Vysoka/phd?/deal/NSE/CMakeLists.txt).

## Run

Run the executable from the `build/` directory so the mesh path `../pipe.msh` resolves correctly:

```bash
cd build
./step57base
```

## Output Directory

You can choose where result files are written with `--output-dir <path>`.
If the directory does not exist, it is created automatically.

Example:

```bash
cd build
./step57base --output-dir ../results/run_01
```

If `--output-dir` is omitted, output files are written to the current working directory.
