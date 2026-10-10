# Local verification record

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
