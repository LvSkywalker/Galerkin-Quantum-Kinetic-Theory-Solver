# Third-party software and data

The root MIT license applies to original project code and documentation. It does not replace the terms of bundled third-party components or external datasets.

## Eigen 3.4.0

The bundled `third_party/eigen` header trees retain upstream copyright and license notices. Most Eigen code is MPL 2.0; components with other terms are identified in their source headers and `COPYING.*` files. The bundled KissFFT backend retains its upstream notice. See [the dependency record](third_party/README.md).

## nlohmann/json 3.11.3

The unmodified single header `third_party/nlohmann/json.hpp` includes its copyright and MIT license notice.

## Hilker et al. dataset

The experimental CSVs, original archive and readme under `references/hilker2022/` come from Hilker et al., *Research data supporting First and Second Sound in a Compressible 3D Bose Fluid*, DOI [10.17863/CAM.83837](https://doi.org/10.17863/CAM.83837), licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The full author list, download provenance and archive checksum are in `ATTRIBUTION.json`.

The application selects five observations, converts units and generates model predictions and adapted figures. These changes are documented in [the application note](docs/APPLICATIONS.md). Experimental data and adapted plots must retain the source attribution; they are not original measurements by this project's author.

## Scientific methods

Uehling--Uhlenbeck kinetic theory, Sonine--Laguerre/Galerkin methods, Gaussian quadrature and the variational benchmark of Wu are established work cited in the README and technical note. This project does not claim their invention. Frozen numerical fixtures from the development prototype are documented under `tests/reference/`.
