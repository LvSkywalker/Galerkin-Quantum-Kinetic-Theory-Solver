#include "qtransport/core.hpp"
#include <cmath>
#include <stdexcept>

namespace qt {
Json shear_wave(const Json &r, double k, const Vec &supplied_times) {
    require(std::isfinite(k) && k > 0, "Wave number must be positive");
    int n = r.at("basis_size");
    Equilibrium eq(r.at("statistics"), r.at("fugacity"), r.at("equilibrium_order"));
    Vec shifts = vector(r.at("thermal_matching_shifts"));
    int size = 2 * n + 1;
    Mat h = Mat::Zero(size, size), m = h, s = h;
    for (int i = 0; i < eq.order; ++i) {
        auto [p, q] = eq.profiles(eq.x[i], n, shifts);
        Vec a(n + 1);
        a[0] = 1;
        a.tail(n) = p;
        double v2 = 2 * eq.x[i];
        h.topLeftCorner(n + 1, n + 1).noalias() += eq.w_measure[i] * v2 / 3 * (a * a.transpose());
        s.topRightCorner(n + 1, n).noalias() +=
            eq.w_measure[i] * v2 * v2 / 15 * (a * q.transpose());
    }
    h.bottomRightCorner(n, n) = matrix(r.at("shear_mass_matrix"));
    s.bottomLeftCorner(n, n + 1) = s.topRightCorner(n + 1, n).transpose();
    m.block(1, 1, n, n) = matrix(r.at("heat_collision_matrix"));
    m.bottomRightCorner(n, n) = matrix(r.at("shear_collision_matrix"));
    Eigen::LLT<Mat> factor(h);
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("Wave mass matrix not positive definite");
    Mat lower = factor.matrixL();
    Mat inverse = lower.triangularView<Eigen::Lower>().solve(Mat::Identity(size, size));
    CMat generator = -inverse.cast<std::complex<double>>() *
                     (m.cast<std::complex<double>>() +
                      std::complex<double>(0, k) * s.cast<std::complex<double>>()) *
                     inverse.transpose().cast<std::complex<double>>();
    Eigen::ComplexEigenSolver<CMat> eig(generator);
    if (eig.info() != Eigen::Success)
        throw std::runtime_error("Wave eigensolver failed");
    CVec initial = lower.row(0).transpose().cast<std::complex<double>>();
    CVec weights = eig.eigenvectors().fullPivLu().solve(initial);
    int slow = 0;
    for (int i = 1; i < size; ++i)
        if (eig.eigenvalues()[i].real() > eig.eigenvalues()[slow].real())
            slow = i;
    double mu = r.at("viscosity");
    Vec times = supplied_times.size() ? supplied_times : linspace(0, 4 / (mu * k * k), 151);
    require(times.size() >= 2 && times[0] == 0, "Time grid must start at zero");
    Vec amplitude(times.size()), entropy(times.size()), ns(times.size());
    double imaginary = 0;
    for (int it = 0; it < times.size(); ++it) {
        require(std::isfinite(times[it]) && (it == 0 || times[it] > times[it - 1]),
                "Time grid must increase");
        CVec y = eig.eigenvectors() *
                 (weights.array() * (times[it] * eig.eigenvalues().array()).exp()).matrix();
        CVec c = inverse.transpose().cast<std::complex<double>>() * y;
        std::complex<double> velocity = (h.row(0).cast<std::complex<double>>() * c)(0);
        amplitude[it] = velocity.real();
        imaginary = std::max(imaginary, std::abs(velocity.imag()));
        entropy[it] = .5 * y.squaredNorm();
        ns[it] = std::exp(-mu * k * k * times[it]);
    }
    auto rate = eig.eigenvalues()[slow];
    double inferred = -rate.real() / (k * k);
    return {
        {"wave_number", k},
        {"time", array(times)},
        {"velocity_amplitude", array(Vec(amplitude / amplitude[0]))},
        {"velocity_imaginary_max_abs", imaginary},
        {"navier_stokes_amplitude", array(ns)},
        {"normalized_quadratic_entropy", array(Vec(entropy / entropy[0]))},
        {"slow_eigenvalue_real", rate.real()},
        {"slow_eigenvalue_imaginary", rate.imag()},
        {"viscosity_from_slow_mode", inferred},
        {"relative_difference_from_transport_viscosity", inferred / mu - 1},
        {"maximum_eigenvalue_real", rate.real()},
        {"note",
         "Internal projected kinetic/hydrodynamic consistency check; not independent validation"}};
}
CVec sound_modes(double qe, double qk, double ratio) {
    require(std::isfinite(qe) && std::isfinite(qk) && qe >= 0 && qk >= 0 && ratio > 0 && ratio <= 1,
            "Invalid sound parameters");
    // Companion matrix of Hilker Eq. S52 for C=omega/(k*c_ad).
    CMat companion = CMat::Zero(3, 3);
    companion(1, 0) = 1;
    companion(2, 1) = 1;
    companion(0, 2) = std::complex<double>(0, qk * ratio);
    companion(1, 2) = 1 + qk * qe;
    companion(2, 2) = std::complex<double>(0, -qk - qe);
    Eigen::ComplexEigenSolver<CMat> eig(companion, false);
    if (eig.info() != Eigen::Success)
        throw std::runtime_error("Sound roots failed");
    return eig.eigenvalues();
}
Json sound_prediction(const Json &s, double eta, double kappa, double length) {
    require(eta >= 0 && kappa > 0 && length > 0, "Invalid transport or length");
    double k = pi / length, rho = s.at("rho"), cv = s.at("cv"), cp = s.at("cp"),
           cad = s.at("sound_speed");
    double qe = 4 * eta * k / (3 * rho * cad), qk = kappa * k / (cv * cad);
    CVec roots = sound_modes(qe, qk, cv / cp);
    int imax = 0;
    Json values = Json::array();
    for (int i = 0; i < 3; ++i) {
        values.push_back({roots[i].real(), roots[i].imag()});
        if (roots[i].real() > roots[imax].real())
            imax = i;
    }
    auto c = roots[imax];
    Json out = {{"q_eta", qe},
                {"q_kappa", qk},
                {"prandtl", cp * eta / (rho * kappa)},
                {"cv_over_cp", cv / cp},
                {"wavenumber", k},
                {"underdamped", c.real() > 1e-9},
                {"roots", values},
                {"damping_stokes_kirchhoff",
                 k * k / 2 * (4 * eta / (3 * rho) + kappa * (1 / cv - 1 / cp))}};
    if (c.real() > 1e-9) {
        out["speed"] = std::abs(c) * cad;
        out["phase_speed"] = c.real() * cad;
        out["damping"] = -c.imag() * k * cad;
        out["oscillation_frequency"] = c.real() * k * cad;
        out["speed_over_hf"] = std::abs(c);
        out["damping_ratio"] = -c.imag() / std::abs(c);
    }
    return out;
}
CMat longitudinal_matrix(const Json &s, double eta, double kappa, double k) {
    double n = s.at("n"), t = s.at("temperature"), rho = s.at("rho"), cv = s.at("cv"),
           c = s.at("sound_speed"), pn = s.at("pressure_density_derivative"),
           pt = s.at("pressure_temperature_derivative");
    CMat m = CMat::Zero(3, 3);
    std::complex<double> I(0, 1);
    m(0, 1) = -I;
    m(1, 0) = -I * n * pn / (rho * c * c);
    m(1, 1) = -4 * eta * k / (3 * rho * c);
    m(1, 2) = -I * t * pt / (rho * c * c);
    m(2, 1) = -I * pt / cv;
    m(2, 2) = -kappa * k / (cv * c);
    return m;
}
} // namespace qt
