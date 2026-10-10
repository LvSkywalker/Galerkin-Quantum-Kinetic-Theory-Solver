#include "qtransport/core.hpp"
#include <cmath>

namespace qt {
namespace {
struct PairFeature {
    Mat a, b;
    double den, energy;
    Eigen::Vector3d momentum;
};
PairFeature pair_feature(const Equilibrium &eq, int n, const Vec &shifts, double center, double q,
                         double cosine, double phi) {
    Eigen::Vector3d relative(q * std::sqrt(std::max(0.0, 1 - cosine * cosine)) * std::cos(phi),
                             q * std::sqrt(std::max(0.0, 1 - cosine * cosine)) * std::sin(phi),
                             q * cosine);
    Eigen::Vector3d v = Eigen::Vector3d(0, 0, center) + relative,
                    v1 = Eigen::Vector3d(0, 0, center) - relative;
    PairFeature out{Mat::Zero(n, 3), Mat::Zero(n, 6), 1, 0, v + v1};
    for (const auto &vel : {v, v1}) {
        double energy = vel.squaredNorm() / 2;
        auto [p, s] = eq.profiles(energy, n, shifts);
        Eigen::Matrix<double, 6, 1> tensor;
        double third = vel.squaredNorm() / 3;
        tensor << vel[0] * vel[0] - third, vel[1] * vel[1] - third, vel[2] * vel[2] - third,
            std::sqrt(2.) * vel[0] * vel[1], std::sqrt(2.) * vel[0] * vel[2],
            std::sqrt(2.) * vel[1] * vel[2];
        out.a.noalias() += p * vel.transpose();
        out.b.noalias() += s * tensor.transpose();
        out.den *= eq.denominator(energy);
        out.energy += energy;
    }
    return out;
}
} // namespace
Brackets brackets(const Equilibrium &eq, const Options &o) {
    require(o.basis >= 1 && o.radial >= 2 && o.angular >= 2 && o.azimuth >= 3,
            "Invalid basis/quadrature size");
    require(std::isfinite(o.gamma) && o.gamma > -3 && std::isfinite(o.strength) && o.strength > 0 &&
                std::isfinite(o.diameter) && o.diameter > 0,
            "Invalid collision kernel");
    auto rx = laguerre_rule(o.radial, .5), ry = laguerre_rule(o.radial, (1 + o.gamma) / 2),
         polar = legendre_rule(o.angular);
    const int n = o.basis;
    Vec shift = eq.shifts(n);
    Mat heat = Mat::Zero(n, n), shear = heat, da(n, 3), db(n, 6);
    double energy_error = 0, momentum_error = 0;
    // Pair features are cached per radius pair. Only the relative orientation
    // and outer products remain in the innermost collision quadrature loop.
    for (int ix = 0; ix < o.radial; ++ix)
        for (int iy = 0; iy < o.radial; ++iy) {
            double center = std::sqrt(rx.x[ix]), q = std::sqrt(ry.x[iy]);
            std::vector<PairFeature> incoming, outgoing;
            incoming.reserve(o.angular);
            outgoing.reserve(o.angular * o.azimuth);
            for (int ia = 0; ia < o.angular; ++ia)
                incoming.push_back(pair_feature(eq, n, shift, center, q, polar.x[ia], 0));
            for (int ib = 0; ib < o.angular; ++ib)
                for (int ip = 0; ip < o.azimuth; ++ip)
                    outgoing.push_back(pair_feature(eq, n, shift, center, q, polar.x[ib],
                                                    2 * pi * ip / o.azimuth));
            double prefactor = 4 * std::pow(2., o.gamma) * o.strength * rx.w[ix] * ry.w[iy] /
                               (4 * o.azimuth * eq.g3_over_z * eq.g3_over_z);
            for (int ia = 0; ia < o.angular; ++ia)
                for (int ib = 0; ib < o.angular; ++ib)
                    for (int ip = 0; ip < o.azimuth; ++ip) {
                        const auto &pre = incoming[ia];
                        const auto &post = outgoing[ib * o.azimuth + ip];
                        da = pre.a - post.a;
                        db = pre.b - post.b;
                        double measure =
                            prefactor * polar.w[ia] * polar.w[ib] / (pre.den * post.den);
                        energy_error = std::max(energy_error, std::abs(pre.energy - post.energy));
                        momentum_error = std::max(
                            momentum_error, (pre.momentum - post.momentum).cwiseAbs().maxCoeff());
                        for (int i = 0; i < n; ++i)
                            for (int j = 0; j <= i; ++j) {
                                heat(i, j) += measure / 3 * da.row(i).dot(da.row(j));
                                shear(i, j) += measure / 10 * db.row(i).dot(db.row(j));
                            }
                    }
        }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j) {
            heat(j, i) = heat(i, j);
            shear(j, i) = shear(i, j);
        }
    return {heat,
            shear,
            {{"collision_energy_abs_error", energy_error},
             {"collision_momentum_abs_error", momentum_error},
             {"collision_quadruples",
              int64_t(o.radial) * o.radial * o.angular * o.angular * o.azimuth}}};
}
Json wu_reference(const std::string &stat, double z, int radial, int angular) {
    Equilibrium eq(stat, z);
    auto polar = legendre_rule(angular), ry = laguerre_rule(radial, 3);
    double integrals[2] = {0, 0};
    for (int sector = 0; sector < 2; ++sector) {
        auto rx = laguerre_rule(radial, sector == 0 ? 1.5 : .5);
        long double total = 0;
        for (int ix = 0; ix < radial; ++ix)
            for (int iy = 0; iy < radial; ++iy)
                for (int ia = 0; ia < angular; ++ia)
                    for (int ib = 0; ib < angular; ++ib) {
                        double a = (polar.x[ia] + 1) / 2, b = (polar.x[ib] + 1) / 2, x = rx.x[ix],
                               y = ry.x[iy], den = 1;
                        for (double cosine : {a, b})
                            for (int sign : {-1, 1})
                                den *= eq.denominator(
                                    (x + y + sign * 2 * std::sqrt(x * y) * cosine) / 2);
                        double angle = sector == 0 ? a * a + b * b - 2 * a * a * b * b
                                                   : 1 + a * a + b * b - 3 * a * a * b * b;
                        total +=
                            rx.w[ix] * ry.w[iy] * polar.w[ia] * polar.w[ib] * angle / (16 * den);
                    }
        integrals[sector] = double(total);
    }
    double g[3];
    for (int i = 0; i < 3; ++i) {
        double nu = 1.5 + i;
        auto rule = laguerre_rule(std::max(96, radial), nu - 1);
        g[i] = 0;
        for (int k = 0; k < rule.x.size(); ++k)
            g[i] += rule.w[k] / eq.denominator(rule.x[k]) / std::tgamma(nu);
    }
    return {
        {"viscosity_ratio", g[1] * g[1] * std::sqrt(pi) / integrals[1]},
        {"conductivity_ratio",
         std::pow(3.5 * g[2] - 2.5 * g[1] * g[1] / g[0], 2) * std::sqrt(pi) / 2 / integrals[0]},
        {"radial_order", radial},
        {"angular_order", angular},
        {"reference", "Wu arXiv:1810.03090v1 Eq. (50), one-profile normalized variational curves"}};
}
} // namespace qt
