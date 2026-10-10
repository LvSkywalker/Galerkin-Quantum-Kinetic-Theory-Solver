#include "qtransport/core.hpp"
#include <cmath>
#include <stdexcept>

namespace qt {
Json relaxation(const Mat &h, const Mat &m, const Vec &d, const Vec &times) {
    Eigen::GeneralizedSelfAdjointEigenSolver<Mat> eig(m, h);
    if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() <= 0)
        throw std::runtime_error("Nonpositive or unresolved projected relaxation spectrum");
    Vec rates = eig.eigenvalues(), amp = eig.eigenvectors().transpose() * d;
    Vec current(times.size()), entropy(times.size());
    for (int i = 0; i < times.size(); ++i) {
        require(std::isfinite(times[i]) && times[i] >= 0, "Invalid time");
        current[i] = (amp.array().square() * (-times[i] * rates.array()).exp()).sum();
        entropy[i] = .5 * (amp.array().square() * (-2 * times[i] * rates.array()).exp()).sum();
    }
    return {{"rates", array(rates)},
            {"time", array(times)},
            {"normalized_current", array(Vec(current / current[0]))},
            {"normalized_quadratic_entropy", array(Vec(entropy / entropy[0]))},
            {"response_time_integral", (amp.array().square() / rates.array()).sum()}};
}
Json solve(const Options &o) {
    Equilibrium eq(o.statistics, o.z, o.equilibrium);
    auto mom = eq.moments(o.basis);
    auto b = brackets(eq, o);
    Eigen::LLT<Mat> lh(b.heat), ls(b.shear);
    if (lh.info() != Eigen::Success || ls.info() != Eigen::Success)
        throw std::runtime_error("Collision bracket is not positive definite");
    Vec ah = lh.solve(mom.dh), bs = ls.solve(mom.ds);
    double conductivity = mom.dh.dot(ah), viscosity = mom.ds.dot(bs);
    double scale = 4 * o.strength * std::pow(2., o.gamma - 1) * std::tgamma((o.gamma + 7) / 2) / 6;
    double mu1 = 5 / (16 * scale * std::sqrt(pi)), kappa1 = 75 / (64 * scale * std::sqrt(pi));
    double moment_k = 0, moment_mu = 0, compat = 0;
    for (int i = 0; i < eq.order; ++i) {
        auto [p, q] = eq.profiles(eq.x[i], o.basis, mom.shifts);
        double v2 = 2 * eq.x[i];
        moment_k += eq.w_measure[i] * p.dot(ah) * v2 * v2 / 6;
        moment_mu += eq.w_measure[i] * q.dot(bs) * v2 * v2 / 15;
        compat += eq.w_measure[i] * v2 / 3 * (eq.x[i] - 2.5 * eq.pressure);
    }
    b.diagnostics["thermal_matching_max_abs"] = mom.matching.cwiseAbs().maxCoeff();
    b.diagnostics["source_compatibility_abs"] = std::abs(compat);
    b.diagnostics["conductivity_moment_relative_difference"] =
        std::abs(moment_k / conductivity - 1);
    b.diagnostics["viscosity_moment_relative_difference"] = std::abs(moment_mu / viscosity - 1);
    b.diagnostics["heat_solve_relative_residual"] = (b.heat * ah - mom.dh).norm() / mom.dh.norm();
    b.diagnostics["shear_solve_relative_residual"] = (b.shear * bs - mom.ds).norm() / mom.ds.norm();
    for (const std::string sector : {"heat", "shear"}) {
        Eigen::SelfAdjointEigenSolver<Mat> eigen(sector == "heat" ? b.heat : b.shear,
                                                 Eigen::EigenvaluesOnly);
        if (eigen.info() != Eigen::Success)
            throw std::runtime_error("Bracket eigenvalue diagnostic failed");
        b.diagnostics[sector + "_bracket_min_eigenvalue"] = eigen.eigenvalues().minCoeff();
        b.diagnostics[sector + "_matrix_condition_number"] =
            eigen.eigenvalues().maxCoeff() / eigen.eigenvalues().minCoeff();
    }
    Vec times = linspace(0, 2 / (4 * o.strength), 101);
    return {{"statistics", o.statistics},
            {"fugacity", o.z},
            {"basis_size", o.basis},
            {"radial_order", o.radial},
            {"angular_order", o.angular},
            {"azimuthal_order", o.azimuth},
            {"equilibrium_order", o.equilibrium},
            {"diameter", o.diameter},
            {"kernel_gamma", o.gamma},
            {"kernel_strength", o.strength},
            {"units", "n=m=k_B=T=1; isotropic B=C*g_rel^gamma"},
            {"pressure", eq.pressure},
            {"degeneracy_per_state", eq.z * eq.g3_over_z},
            {"occupation_conversion_A", eq.occupation_conversion},
            {"viscosity", viscosity},
            {"conductivity", conductivity},
            {"viscosity_over_classical_first_sonine", viscosity / mu1},
            {"conductivity_over_classical_first_sonine", conductivity / kappa1},
            {"classical_first_sonine_viscosity", mu1},
            {"classical_first_sonine_conductivity", kappa1},
            {"thermal_coefficients", array(ah)},
            {"viscous_coefficients", array(bs)},
            {"thermal_matching_shifts", array(mom.shifts)},
            {"heat_collision_matrix", array(b.heat)},
            {"shear_collision_matrix", array(b.shear)},
            {"heat_mass_matrix", array(mom.hh)},
            {"shear_mass_matrix", array(mom.hs)},
            {"diagnostics", b.diagnostics},
            {"heat_relaxation", relaxation(mom.hh, b.heat, mom.dh, times)},
            {"shear_relaxation", relaxation(mom.hs, b.shear, mom.ds, times)},
            {"implementation", "C++17 / Eigen"}};
}
} // namespace qt
