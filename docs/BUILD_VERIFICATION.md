# Local verification record

## Version 0.4.1 - 10 October 2026

Fresh Release build with Apple Clang 21.0.0 and CMake on macOS arm64:

| Check | Result |
|---|---|
| CMake build | Passed, including the C++ report exporter |
| CTest | 3/3 passed in 2.13 s; 24 physics/regression groups and 2 CLI checks |
| Production reproduction | 36 transport calculations, five N=1 references, basis/quadrature refinement and linear dynamics; 26.46 s |
| Numerical preservation | Transport, reference, convergence, dynamics and saved Hilker arrays are byte-identical to the downloaded 0.4.0 publication |
| Physical-core preservation | Collision, equilibrium, transport, thermodynamics and dynamics source files are unchanged |
| Frozen-reference audit | Passed: 36 transport cases and five saved sound predictions |
| Report-export verification | Independently recomputed all four standardized RMS values from the selected observations |
| Figures and note | Rebuilt both revised figures and the 16-page note; rendered and visually reviewed |
| README/document links | All checked relative targets exist |

The Fourier diagnostic and its two associated tests were removed; no transport or sound test was removed. The Hilker collision calculation was not rerun: the existing saved predictions are unchanged, and this revision adds their explicit classical-transport comparison. The updated exporter writes per-observation residuals to `docs/technical_note/data/hilker_metrics.json`.

These results are local verification of 0.4.1. GitHub Actions must run again after uploading this revision. The earlier public 0.4.0 workflow completed successfully, but that does not certify an unpublished revision. Source-archive SHA-256 hashes are generated only inside release ZIPs.

## Previous verification: version 0.4.0


Verified on 9 October 2026 for the prepared 0.4.0 source tree.

Environment: macOS/Darwin arm64, Apple Clang 21.0.0 (clang-2100.1.1.101), CMake 3.31. Both builds used the macOS Command Line Tools SDK. This is local evidence; Linux and remote GitHub Actions results are not asserted here.

| Check | Result |
|---|---|
| Fresh CMake Release build | Completed |
| Release CTest | 3/3 passed: 26 physics/regression groups plus 2 CLI checks |
| Debug build with ASan + UBSan, `-O1 -g`, frame pointers | Completed |
| Sanitized CTest | 3/3 passed, no sanitizer findings |
| Production transport reproduction | 36 calculations; five N=1 benchmarks; basis and quadrature refinement |
| Optional Hilker production reproduction | Five selected observations and recorded sensitivity calculations |
| Frozen-reference audit | Passed for 36 transport cases and five sound predictions |
| Technical-note build | Tectonic, 15 pages, no TeX warnings in the final build |
| Document review | Rendered all pages; checked layout, figures and references |

The Release test suite took 2.13 s locally; the sanitized suite took 44.87 s. These are test-suite timings on this execution, not portable performance benchmarks or speedup claims. An initial unoptimized sanitizer run was stopped and rerun with the documented `-O1` sanitizer configuration; the results above refer to the completed run.

The largest one-profile discrepancy was 0.00302269%. The largest displayed N=4 to N=5 change was 0.313386%, normalized to N=4. Frozen-reference relative differences were at most 7.21e-12 for transport and 2.65e-11 for sound predictions. Their respective meanings are explained in [VALIDATION.md](VALIDATION.md).

The fresh structured evidence is in `results/cpp/test_report.json`, `manifest.json`, `benchmark.json`, `convergence.json`, `regression_audit.json` and `hilker/manifest.json`. Source-archive SHA-256 hashes are recorded in `PACKAGE_MANIFEST.json` when packaging.
