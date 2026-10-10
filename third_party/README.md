# Bundled dependencies

Dependencies are pinned and bundled so a build does not require network access.

- **Eigen 3.4.0**, downloaded from https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz. The `Eigen` and `unsupported/Eigen` header trees are unchanged. Upstream license notices are in `eigen/COPYING.*`. Dense algebra, symmetric/generalized/complex eigensolvers and the bundled KissFFT backend are used. Unused upstream documentation, examples and build scripts are omitted.
- **nlohmann/json 3.11.3**, unchanged single header from https://raw.githubusercontent.com/nlohmann/json/v3.11.3/single_include/nlohmann/json.hpp. The MIT license and copyright notice are included in the header.

The project implements the quadrature construction, equilibrium, collision integration, physical models, discretization and validation. It does not claim these third-party algebra/FFT/JSON implementations as original work.
