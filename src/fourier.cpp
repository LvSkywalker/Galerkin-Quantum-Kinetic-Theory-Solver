#include "qtransport/core.hpp"
#include <cmath>
#include <unsupported/Eigen/FFT>

namespace qt {
namespace {
using Field = std::vector<std::complex<double>>;
size_t index(int i, int j, int k, int n) {
    return (size_t(i) * n + j) * n + k;
}
Field transform(Field values, int n, bool inverse) {
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> line(n), out(n);
    for (int axis = 0; axis < 3; ++axis)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                auto at = [&](int k) {
                    return axis == 0 ? index(k, i, j, n)
                                     : (axis == 1 ? index(i, k, j, n) : index(i, j, k, n));
                };
                for (int k = 0; k < n; ++k)
                    line[k] = values[at(k)];
                if (inverse)
                    fft.inv(out, line);
                else
                    fft.fwd(out, line);
                for (int k = 0; k < n; ++k)
                    values[at(k)] = out[k];
            }
    return values;
}
std::pair<Field, Field> sample(const Equilibrium &eq, const Vec &b, double width, int n,
                               std::array<double, 3> offset) {
    size_t size = size_t(n) * n * n;
    Field f(size), g(size);
    double dx = 2 * width / n;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                double x = dx * (i - n / 2 + offset[0]), y = dx * (j - n / 2 + offset[1]),
                       z = dx * (k - n / 2 + offset[2]), e = (x * x + y * y + z * z) / 2,
                       den = eq.denominator(e);
                size_t a = index(i, j, k, n);
                f[a] = std::exp(-e) / (eq.normalization * den);
                g[a] = f[a] / den * laguerre_values(e, int(b.size()), 2.5).dot(b) * x * y;
            }
    return {f, g};
}
Field coefficients(const Field &values, int n) {
    Field shifted(values.size());
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                shifted[index(i, j, k, n)] =
                    values[index((i + n / 2) % n, (j + n / 2) % n, (k + n / 2) % n, n)];
    auto out = transform(std::move(shifted), n, false);
    for (auto &v : out)
        v /= double(values.size());
    return out;
}
Field reconstruct(Field coeff, int n, std::array<double, 3> offset) {
    auto freq = [&](int k) { return k <= n / 2 ? k : k - n; };
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                coeff[index(i, j, k, n)] *= std::exp(std::complex<double>(
                    0, 2 * pi * (freq(i) * offset[0] + freq(j) * offset[1] + freq(k) * offset[2]) /
                           n));
    auto raw = transform(std::move(coeff), n, true);
    Field out(raw.size());
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                out[index(i, j, k, n)] =
                    double(raw.size()) * raw[index((i + (n + 1) / 2) % n, (j + (n + 1) / 2) % n,
                                                   (k + (n + 1) / 2) % n, n)];
    return out;
}
} // namespace
Json representation(const Json &r, double width, int n) {
    require(std::isfinite(width) && width > 0 && n >= 3 && n % 2 == 1,
            "Fourier grid requires positive width and odd points >=3");
    Equilibrium eq(r.at("statistics"), r.at("fugacity"), r.at("equilibrium_order"));
    Vec b = vector(r.at("viscous_coefficients"));
    auto [f, g] = sample(eq, b, width, n, {0, 0, 0});
    auto target = sample(eq, b, width, n, {.37, .23, .41});
    double dx = 2 * width / n, volume = dx * dx * dx, density = 0, pressure = 0, stress = 0;
    std::array<double, 5> matching{0, 0, 0, 0, 0};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                double x = dx * (i - n / 2), y = dx * (j - n / 2), z = dx * (k - n / 2),
                       v2 = x * x + y * y + z * z;
                auto a = index(i, j, k, n);
                density += f[a].real() * volume;
                pressure += v2 * f[a].real() * volume / 3;
                stress += x * y * g[a].real() * volume;
                std::array<double, 5> invariants{1, x, y, z, v2};
                for (int l = 0; l < 5; ++l)
                    matching[l] += invariants[l] * g[a].real() * volume;
            }
    double matching_error = 0;
    for (double x : matching)
        matching_error = std::max(matching_error, std::abs(x));
    Json out = {{"statistics", eq.statistics},
                {"fugacity", eq.z},
                {"basis_size", b.size()},
                {"half_width", width},
                {"points_per_axis", n},
                {"velocity_spacing", dx},
                {"density_relative_error", density - 1},
                {"pressure_relative_error", pressure / eq.pressure - 1},
                {"shear_moment_relative_error", stress / r.at("viscosity").get<double>() - 1},
                {"reference_viscosity_galerkin", r.at("viscosity")},
                {"shear_matching_max_abs", matching_error},
                {"note", "Velocity Fourier representation only; no Fourier collision solver"}};
    for (int sector = 0; sector < 2; ++sector) {
        const auto &values = sector ? g : f;
        const auto &exact = sector ? target.second : target.first;
        std::string label = sector ? "shear" : "equilibrium";
        auto coef = coefficients(values, n), back = reconstruct(coef, n, {0, 0, 0}),
             off = reconstruct(coef, n, {.37, .23, .41});
        double norm = 0, shift_norm = 0, err = 0, shift_err = 0, peak = 0, edge = 0, imag = 0,
               min = off[0].real();
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    auto a = index(i, j, k, n);
                    norm += std::norm(values[a]);
                    shift_norm += std::norm(exact[a]);
                    err += std::norm(back[a] - values[a]);
                    shift_err += std::norm(off[a] - exact[a]);
                    peak = std::max(peak, std::abs(values[a]));
                    imag = std::max(imag, std::abs(off[a].imag()));
                    min = std::min(min, off[a].real());
                    if (i == 0 || j == 0 || k == 0 || i == n - 1 || j == n - 1 || k == n - 1)
                        edge = std::max(edge, std::abs(values[a]));
                }
        out[label + "_fft_roundtrip_relative_l2"] = std::sqrt(err / norm);
        out[label + "_offgrid_relative_l2"] = std::sqrt(shift_err / shift_norm);
        out[label + "_boundary_layer_relative_max"] = edge / peak;
        out[label + "_imaginary_relative_max"] = imag / peak;
        if (!sector)
            out["equilibrium_interpolated_min"] = min;
    }
    return out;
}
} // namespace qt
