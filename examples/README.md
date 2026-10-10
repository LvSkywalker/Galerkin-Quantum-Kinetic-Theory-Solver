# Small runnable examples

Build first and run these commands from the repository root. Each transport command creates a JSON file under the ignored `runs/` directory, leaving the saved production scan in place.

## Classical first-profile benchmark

```sh
build/qtransport transport --statistics classical --basis 1 --output runs/classical-n1.json
```

Expected dimensionless values: viscosity approximately 0.1763092449 and conductivity approximately 0.6611596682.

## Bose transport with four profiles

```sh
build/qtransport transport --statistics bose --fugacity 0.5 --basis 4 \
  --radial 32 --angular 24 --azimuth 24 --output runs/bose-z0.5-n4.json
```

Compare to the classical calculation at the same N and integration settings when forming a quantum/classical ratio.

## Separate one-profile reference integrals

```sh
build/qtransport transport --statistics bose --fugacity 0.5 --basis 1 \
  --radial 32 --angular 24 --azimuth 24 --output runs/bose-z0.5-n1.json
build/qtransport reference --statistics bose --fugacity 0.5 \
  --radial 64 --angular 40 --output runs/wu-bose-z0.5.json
```

Compare `viscosity_over_classical_first_sonine` with the reference's `viscosity_ratio`, and similarly for conductivity. This is the same N=1 approximation evaluated through different integrals.

## Change the collision frequency scale

```sh
build/qtransport transport --statistics fermi --fugacity 1 --basis 4 \
  --gamma 1 --strength 0.25 --output runs/fermi-c025.json
build/qtransport transport --statistics fermi --fugacity 1 --basis 4 \
  --gamma 1 --strength 0.5 --output runs/fermi-c050.json
```

Doubling the kernel strength halves the transport coefficients and doubles the homogeneous relaxation rates, at otherwise identical settings. This scaling is also checked by the automated suite.

## Follow a transverse mode

```sh
build/qtransport wave --input runs/bose-z0.5-n4.json --k 0.1 \
  --output runs/bose-shear-k01.json
```

This projects streaming and collisions onto a finite velocity basis. Its small-k limit is checked against steady viscosity. Homogeneous relaxation already appears in every `transport` result under `heat_relaxation` and `shear_relaxation`.
