# Changelog

## 0.4.0 - 2026-10-09

- Prepared a standalone C++ repository with build instructions, examples, citation metadata and dependency/data attribution.
- Made the technical note independently readable, including the kinetic equation, units, constitutive moments and the exact meaning of the one-profile benchmark.
- Moved prototype regression fixtures under `tests/reference/` and documented their provenance.
- Added CLI help and CTest smoke checks; configured sanitizer builds with debug information and light optimization.
- Reproduced the transport and optional sound calculations after reorganization. The physical model, collision integrals and Galerkin algorithms are unchanged from the preceding C++ implementation.

## Development provenance

The first implementation was a Python prototype. The subsequent C++17 implementation retained numerical fixtures to check matrices, coefficients and transients during migration. Those fixtures remain part of the regression suite; Python is not required to build or execute this repository.
