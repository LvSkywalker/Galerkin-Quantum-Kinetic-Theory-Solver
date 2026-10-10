# Reproducing the calculations

Commands below run from the repository root. A compiler with C++17 support is required. Dependencies are bundled; no package download is part of the solver build.

## Build and tests

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

The compiler-only fallback is `sh build.sh`, followed by `build/qtransport test --root .`. On macOS it selects the separate Command Line Tools compiler when available; `CXX` can select another C++17 compiler.

For a local sanitizer build with Clang or GCC:

```sh
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug -DQTRANSPORT_SANITIZERS=ON
cmake --build build-sanitize --parallel 2
ctest --test-dir build-sanitize --output-on-failure
```

The GitHub workflow specifies Linux/macOS Debug and Release jobs and a Linux sanitizer job. A configured workflow is not evidence of a successful remote run; its status is established when the repository is pushed.
The sanitizer option adds `-O1 -g` and frame pointers alongside AddressSanitizer and UndefinedBehaviorSanitizer, keeping the quadrature-heavy checks practical while retaining debugging information.

## One calculation

```sh
build/qtransport transport --statistics bose --fugacity 0.5 --basis 4 \
  --radial 32 --angular 24 --azimuth 24 --output runs/bose.json
```

`--basis` sets the number of radial profiles. The integration orders independently control radial, polar and azimuthal resolution. The JSON stores all these settings with the numerical results. It contains steady coefficients, collision and Gram matrices, matching shifts, condition numbers, conservation diagnostics and linear relaxation curves. Default settings are convenient for individual runs; use the production workflow below to reproduce the documented evidence.

Units and normalization are defined in [METHODS.md](METHODS.md). In particular, fields named `*_over_classical_first_sonine` use the classical N=1 denominator. For a quantum/classical comparison at the same N, use the ratios in `results/cpp/transport_same_order.csv` or compute both states at that N.

## Published project figures

```sh
build/qtransport reproduce --root .
build/qtransport figures --root .
```

`reproduce` performs 36 transport calculations, the five-state N=1 reference check, basis and quadrature refinement, homogeneous dynamics and projected shear modes. `figures` reads the resulting files and writes SVGs plus `results/cpp/index.html`. It does not rerun the collision integrals. Run scope and elapsed time are stored in `manifest.json`; timing depends on hardware and compiler.

The `--quick` option uses coarser integration and overwrites the corresponding output files. Use it for development only, then rerun the production command before quoting benchmark numbers or committing results.

The optional sound application and regression audit are:

```sh
build/qtransport hilker --root .
build/qtransport figures --root .
build/qtransport audit --root .
```

`audit` compares full transport and Hilker outputs with frozen prototype results in `tests/reference/`. The fixtures' provenance is documented there. Regression agreement is distinct from physical validation.

## Technical note

The [standalone note](technical_note/main.pdf) includes ready-made PDF figures, so compiling the text does not require rebuilding numerical results or installing plotting packages:

```sh
cd docs/technical_note
latexmk -pdf main.tex
```

Tectonic can be used instead: `tectonic main.tex`. LaTeX is only needed for the document, not for the C++ solver.

To update tables and plot data after a production run, from the repository root:

```sh
build/qtransport_report_export . docs/technical_note
```

The exporter requires the transport and Hilker outputs. Rebuild the five figures from their PGFPlots sources, then recompile the note:

```sh
cd docs/technical_note/figures
for f in transport convergence relaxation shear_wave hilker; do
  tectonic "$f.tex"
done
cd ..
tectonic main.tex
```

## Making a source archive

```sh
cmake -P tools/package.cmake
```

The archive excludes build trees, version-control metadata and local `runs/`. A manifest records SHA-256 hashes. Dependency licenses, reference-data attribution, source code and documented results are included.
