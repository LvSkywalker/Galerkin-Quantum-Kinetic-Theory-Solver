# Bundled dependencies

Both dependencies are header-only, pinned and bundled, so the build needs no network access.

- **Eigen 3.4.0** (https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz). The `Eigen` header tree is unchanged; documentation, examples, build scripts and the `unsupported` modules are omitted. License notices are in `eigen/COPYING.*`. Used for dense linear algebra: Cholesky factorizations and symmetric, generalized and complex eigensolvers.
- **nlohmann/json 3.11.3** (https://github.com/nlohmann/json, single header `v3.11.3/single_include/nlohmann/json.hpp`), unchanged. The MIT license notice is in the header. Used for input/output files.
