#include "qtransport/core.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace qt {
namespace {
void check(bool yes, const std::string &msg) {
    if (!yes)
        throw std::runtime_error(msg);
}
void near(double value, double reference, double rtol = 1e-10, double atol = 1e-12) {
    if (!std::isfinite(value) || std::abs(value - reference) > atol + rtol * std::abs(reference))
        throw std::runtime_error("Mismatch: " + std::to_string(value) + " vs " +
                                 std::to_string(reference));
}
void rejects(const std::function<void()> &f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    check(rejected, "Invalid input accepted");
}
void decreasing(const Json &values) {
    double prev = values.at(0);
    near(prev, 1);
    for (const auto &val : values) {
        double v = val;
        check(v >= 0 && v <= prev + 1e-12, "Nonmonotone decay");
        prev = v;
    }
}
Options small() {
    Options o;
    o.radial = 12;
    o.angular = 8;
    o.azimuth = 10;
    return o;
}
} // namespace
int run_tests(const std::string &root) {
    Json report = Json::array();
    int failures = 0;
    auto test = [&](std::string name, const std::function<void()> &body) {
        try {
            body();
            std::cout << "PASS " << name << std::endl;
            report.push_back({{"name", name}, {"passed", true}});
        } catch (const std::exception &e) {
            ++failures;
            std::cout << "FAIL " << name << ": " << e.what() << std::endl;
            report.push_back({{"name", name}, {"passed", false}, {"error", e.what()}});
        }
    };
    test("Gaussian Laguerre moments", [] {
        for (double alpha : {.5, 1., 1.5, 2.5, 3.}) {
            auto q = laguerre_rule(12, alpha);
            for (int k = 0; k < 8; ++k)
                near((q.w.array() * q.x.array().pow(k)).sum(), std::tgamma(alpha + k + 1), 2e-12);
        }
    });
    test("Gaussian Legendre moments", [] {
        auto q = legendre_rule(10);
        for (int k = 0; k < 18; ++k)
            near((q.w.array() * q.x.array().pow(k)).sum(), k % 2 ? 0 : 2. / (k + 1));
    });
    test("Laguerre recurrence", [] {
        auto p = laguerre_values(.7, 3, 1.5);
        near(p[0], 1);
        near(p[1], 1.8);
        near(p[2], (.49 - 7 * .7 + 8.75) / 2);
    });
    test("Equilibrium normalization and matching", [] {
        for (auto pair : std::vector<std::pair<std::string, double>>{
                 {"classical", .1}, {"bose", .5}, {"fermi", 10}}) {
            Equilibrium eq(pair.first, pair.second);
            near(eq.f_measure.sum(), 1);
            check(eq.moments(4).matching.cwiseAbs().maxCoeff() < 1e-10, "Matching violation");
            double compat = 0;
            for (int i = 0; i < eq.order; ++i)
                compat += eq.w_measure[i] * 2 * eq.x[i] * (eq.x[i] - 2.5 * eq.pressure);
            check(std::abs(compat) < 1e-8, "Source incompatible");
        }
    });
    test("Invalid equilibrium inputs", [] {
        for (double z : {0., -1., 1., 1.1})
            rejects([&] { Equilibrium eq("bose", z); });
        rejects([] { Equilibrium eq("unknown", 1); });
    });
    test("Detailed balance and Pauli bound", [] {
        Eigen::Vector3d center(.4, -.1, .8), q(.5, .2, .3),
            qp = q.norm() * Eigen::Vector3d(.6, 0, .8);
        for (std::string stat : {"classical", "bose", "fermi"}) {
            Equilibrium eq(stat, stat == "bose" ? .7 : 10);
            double f[4], r[4];
            std::array<Eigen::Vector3d, 4> velocities{center + q, center - q, center + qp,
                                                      center - qp};
            for (int i = 0; i < 4; ++i) {
                f[i] = eq.density(velocities[i]);
                r[i] = 1 + eq.eps * eq.occupation_conversion * f[i];
                if (stat == "fermi")
                    check(eq.occupation_conversion * f[i] < 1, "Pauli bound");
            }
            near(f[0] * f[1] * r[2] * r[3], f[2] * f[3] * r[0] * r[1], 1e-12, 1e-16);
        }
    });
    test("Classical absolute first-Sonine prefactors", [] {
        auto o = small();
        o.basis = 1;
        auto r = solve(o);
        near(r["viscosity"], 5 / (16 * std::sqrt(pi)));
        near(r["conductivity"], 75 / (64 * std::sqrt(pi)));
        check(r["diagnostics"]["collision_energy_abs_error"].get<double>() < 1e-11,
              "Collision energy");
    });
    test("Power-law classical prefactors", [] {
        for (double gamma : {-2., -1., 0., 1., 2.}) {
            auto o = small();
            o.basis = 1;
            o.gamma = gamma;
            o.strength = .7;
            auto r = solve(o);
            near(r["viscosity_over_classical_first_sonine"], 1);
            near(r["conductivity_over_classical_first_sonine"], 1);
        }
    });
    test("Second-Sonine exact ratios", [] {
        for (double gamma : {-1., 1.}) {
            auto o = small();
            o.gamma = gamma;
            auto r = solve(o);
            near(r["viscosity_over_classical_first_sonine"], 1 + 3. / (gamma == 1 ? 202 : 190));
            near(r["conductivity_over_classical_first_sonine"], 1 + 1. / (gamma == 1 ? 44 : 40));
        }
    });
    test("Kernel strength scaling", [] {
        auto o = small();
        o.statistics = "bose";
        o.z = .5;
        auto a = solve(o);
        o.strength *= 2;
        auto b = solve(o);
        for (std::string key : {"viscosity", "conductivity"})
            near(b[key].get<double>() / a[key].get<double>(), .5);
    });
    test("Invalid collision parameters", [] {
        auto o = small();
        o.gamma = -3;
        rejects([&] { solve(o); });
        o.gamma = 1;
        o.strength = 0;
        rejects([&] { solve(o); });
        o.strength = .25;
        o.basis = 0;
        rejects([&] { solve(o); });
    });
    test("Classical N=4 frozen prototype regression", [] {
        auto o = small();
        o.basis = 4;
        o.angular = 12;
        o.azimuth = 16;
        auto r = solve(o);
        near(r["viscosity"], .17913123631242575);
        near(r["conductivity"], .6777775569049495);
    });
    test("Quantum one-profile Wu independent quadrature", [] {
        for (auto p : std::vector<std::pair<std::string, double>>{{"bose", .5}, {"fermi", 1}}) {
            auto o = small();
            o.statistics = p.first;
            o.z = p.second;
            o.basis = 1;
            o.radial = 16;
            o.angular = 12;
            o.azimuth = 16;
            auto r = solve(o), ref = wu_reference(o.statistics, o.z);
            for (std::string key : {"viscosity", "conductivity"})
                near(r[key + "_over_classical_first_sonine"], ref[key + "_ratio"], 2e-5);
        }
    });
    test("Nested Galerkin convergence", [] {
        double prev_mu = 0, prev_k = 0;
        for (int n : {1, 2, 3}) {
            auto o = small();
            o.basis = n;
            auto r = solve(o);
            double mu = r["viscosity"], k = r["conductivity"];
            check(mu > prev_mu && k > prev_k, "Nested variational monotonicity");
            prev_mu = mu;
            prev_k = k;
        }
    });
    test("Homogeneous entropy and response integral", [] {
        auto o = small();
        o.basis = 4;
        auto r = solve(o);
        for (std::string sector : {"heat", "shear"}) {
            decreasing(r[sector + "_relaxation"]["normalized_quadratic_entropy"]);
            near(r[sector + "_relaxation"]["response_time_integral"],
                 r[sector == "heat" ? "conductivity" : "viscosity"]);
        }
        check(r["diagnostics"]["conductivity_moment_relative_difference"].get<double>() < 1e-10,
              "Moment mismatch");
    });
    test("Shear wave hydrodynamic limit", [] {
        auto o = small();
        o.statistics = "fermi";
        o.z = 1;
        auto r = solve(o), a = shear_wave(r, .05), b = shear_wave(r, .2);
        double err = std::abs(a["relative_difference_from_transport_viscosity"].get<double>());
        check(err < 1e-3 &&
                  err < std::abs(b["relative_difference_from_transport_viscosity"].get<double>()),
              "Hydrodynamic limit");
        check(a["maximum_eigenvalue_real"].get<double>() < 0, "Unstable wave");
        decreasing(a["normalized_quadratic_entropy"]);
    });
    test("Fourier off-grid representation", [] {
        auto o = small();
        o.basis = 2;
        auto r = solve(o), a = representation(r, 7, 17), b = representation(r, 7, 33);
        for (std::string sector : {"equilibrium", "shear"}) {
            check(b[sector + "_fft_roundtrip_relative_l2"].get<double>() < 1e-12, "FFT round trip");
            check(b[sector + "_offgrid_relative_l2"].get<double>() <
                      a[sector + "_offgrid_relative_l2"].get<double>(),
                  "Offgrid refinement");
        }
        check(std::abs(b["density_relative_error"].get<double>()) < 1e-8, "Fourier density");
    });
    test("Fourier input validation", [] {
        auto o = small();
        auto r = solve(o);
        rejects([&] { representation(r, 7, 16); });
        rejects([&] { representation(r, -1, 17); });
    });
    test("Bose series versus Gaussian quadrature", [] {
        for (double z : {.1, .5, .9}) {
            Equilibrium eq("bose", z, 256);
            near(bose_functions(z)[1], z * eq.g3_over_z, 1e-8);
        }
    });
    test("Hartree-Fock thermodynamic derivatives", [] {
        double n = 1e19, t = 200e-9, a = 500 * a0;
        auto s = normal_state(n, t, a), np = normal_state(n * 1.00001, t, a),
             nm = normal_state(n * .99999, t, a), tp = normal_state(n, t * 1.00001, a),
             tm = normal_state(n, t * .99999, a);
        near((np["pressure"].get<double>() - nm["pressure"].get<double>()) / (2e-5 * n),
             s["pressure_density_derivative"], 1e-7, 0);
        near((tp["energy_density"].get<double>() - tm["energy_density"].get<double>()) / (2e-5 * t),
             s["cv"], 1e-7, 0);
    });
    test("Absolute Bose SI prefactors", [] {
        auto o = small();
        o.basis = 1;
        auto r = solve(o);
        double a = 500 * a0, t = 200e-9;
        auto [eta, kappa] = dimensional_transport(r, t, a);
        near(eta, 5 * std::sqrt(mass_k39 * kb * t / pi) / (128 * a * a), 1e-10, 0);
        near(kappa, 75 * kb * std::sqrt(kb * t / (pi * mass_k39)) / (512 * a * a), 1e-10, 0);
    });
    test("Sound cubic versus independent balance matrix", [] {
        auto s = normal_state(1e19, 200e-9, 500 * a0);
        double eta = 1e-14, kappa = 1e-10, l = 50e-6;
        auto p = sound_prediction(s, eta, kappa, l);
        Eigen::ComplexEigenSolver<CMat> eig(longitudinal_matrix(s, eta, kappa, pi / l), false);
        for (const auto &root : p["roots"]) {
            std::complex<double> c(root[0].get<double>(), root[1].get<double>()),
                v = std::complex<double>(0, -1) * c;
            double distance = 1e9;
            for (int i = 0; i < 3; ++i)
                distance = std::min(distance, std::abs(v - eig.eigenvalues()[i]));
            check(distance < 1e-9, "Cubic/matrix mismatch");
        }
    });
    test("Sound small-k Stokes-Kirchhoff limit", [] {
        auto s = normal_state(1e19, 200e-9, 500 * a0);
        auto p = sound_prediction(s, 1e-14, 1e-10, .05);
        near(p["damping"], p["damping_stokes_kirchhoff"], 1e-5);
        near(p["speed"], s["sound_speed"], 1e-5);
    });
    test("Published data units and selected subset", [&] {
        auto rows = observations(root + "/references/hilker2022/Figure_4bc_exp.csv");
        check(rows.size() == 16, "Expected 16 observations");
        int count = 0;
        for (const auto &row : rows) {
            if (row["hydro_selected"].get<bool>())
                ++count;
            if (row["id"] != 13)
                near(row["damping_ratio"],
                     row["damping"].get<double>() /
                         (row["speed"].get<double>() * pi / row["length"].get<double>()));
        }
        check(count == 5, "Expected five selected points");
    });
    test("Hilker frozen state and sound predictions", [&] {
        auto legacy = read_json(root + "/tests/reference/hilker_comparison.json");
        for (const auto &row : legacy)
            if (row.contains("prediction")) {
                auto state = state_from_observation(row["observation"]);
                near(state["fugacity"], row["state"]["fugacity"]);
                auto p = row["prediction"];
                auto actual = sound_prediction(state, p["viscosity"], p["conductivity"],
                                               row["observation"]["length"]);
                near(actual["speed"], p["speed"]);
                near(actual["damping"], p["damping"]);
            }
    });
    test("Full matrices and relaxation versus prototype snapshots", [&] {
        auto legacy = read_json(root + "/tests/reference/dynamics.json");
        for (const auto &row : legacy) {
            auto old = row["transport"];
            Options o;
            o.statistics = old["statistics"];
            o.z = old["fugacity"];
            o.basis = old["basis_size"];
            o.radial = old["radial_order"];
            o.angular = old["angular_order"];
            o.azimuth = old["azimuthal_order"];
            o.equilibrium = old["equilibrium_order"];
            auto r = solve(o);
            for (std::string key : {"viscosity", "conductivity"})
                near(r[key], old[key], 2e-9);
            for (std::string key : {"heat_collision_matrix", "shear_collision_matrix",
                                    "heat_mass_matrix", "shear_mass_matrix"}) {
                Mat a = matrix(r[key]), b = matrix(old[key]);
                check((a - b).norm() / b.norm() < 2e-9, "Matrix regression discrepancy");
            }
            for (std::string sector : {"heat", "shear"})
                for (std::string key :
                     {"rates", "normalized_current", "normalized_quadratic_entropy"}) {
                    Vec a = vector(r[sector + "_relaxation"][key]),
                        b = vector(old[sector + "_relaxation"][key]);
                    check((a - b).norm() / b.norm() < 2e-9, "Relaxation regression discrepancy");
                }
            auto wave = shear_wave(r, .2, vector(row["wave"]["time"]));
            near(wave["viscosity_from_slow_mode"], row["wave"]["viscosity_from_slow_mode"], 2e-8);
            check((vector(wave["velocity_amplitude"]) - vector(row["wave"]["velocity_amplitude"]))
                          .norm() < 1e-7,
                  "Wave regression discrepancy");
        }
    });
    write_json(root + "/results/cpp/test_report.json", {{"implementation", "C++17"},
                                                        {"passed", report.size() - failures},
                                                        {"failed", failures},
                                                        {"tests", report}});
    std::cout << report.size() - failures << " passed, " << failures << " failed\n";
    return failures ? 1 : 0;
}
} // namespace qt
