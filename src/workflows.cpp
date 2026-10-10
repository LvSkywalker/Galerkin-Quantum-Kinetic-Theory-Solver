#include "qtransport/core.hpp"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

namespace qt {
namespace {
Json truncate(const Json &r, int n) {
    Equilibrium eq(r.at("statistics"), r.at("fugacity"), r.at("equilibrium_order"));
    auto m = eq.moments(r.at("basis_size"));
    Json out = r;
    out["basis_size"] = n;
    Mat mh = matrix(r.at("heat_collision_matrix")).topLeftCorner(n, n),
        ms = matrix(r.at("shear_collision_matrix")).topLeftCorner(n, n);
    out["conductivity"] = m.dh.head(n).dot(mh.llt().solve(m.dh.head(n)));
    out["viscosity"] = m.ds.head(n).dot(ms.llt().solve(m.ds.head(n)));
    return out;
}
void csv(const std::string &path, const Json &rows, const std::vector<std::string> &keys) {
    std::ofstream f(path);
    if (!f)
        throw std::runtime_error("Cannot write " + path);
    f << std::setprecision(17);
    for (size_t k = 0; k < keys.size(); ++k)
        f << (k ? "," : "") << keys[k];
    f << '\n';
    for (const auto &r : rows) {
        for (size_t k = 0; k < keys.size(); ++k) {
            if (k)
                f << ',';
            const auto &v = r.at(keys[k]);
            if (v.is_string())
                f << v.get<std::string>();
            else
                f << v;
        }
        f << '\n';
    }
}
Json prediction(const Json &row, const Json &r, bool mf = true) {
    auto state = state_from_observation(row, mf);
    auto [eta, kappa] =
        dimensional_transport(r, row.at("temperature"), row.at("scattering_length"));
    auto out = sound_prediction(state, eta, kappa, row.at("length"));
    out["state"] = state;
    out["viscosity"] = eta;
    out["conductivity"] = kappa;
    return out;
}
} // namespace
Json reproduce(const std::string &root, bool quick) {
    auto start = std::chrono::steady_clock::now();
    std::string out = root + "/results/cpp";
    std::filesystem::create_directories(out);
    Options o;
    o.equilibrium = 128;
    o.radial = quick ? 12 : 24;
    o.angular = quick ? 8 : 20;
    o.azimuth = quick ? 10 : 24;
    Json records = Json::array(), dynamics = Json::array(), benchmarks = Json::array(),
         convergence = Json::array();
    std::map<std::string, Json> full;
    auto key = [](std::string stat, double z, int n) {
        return stat + "_" + std::to_string(z) + "_" + std::to_string(n);
    };
    std::vector<std::pair<std::string, double>> cases = {{"classical", .1}};
    for (double z :
         (quick ? std::vector<double>{.5} : std::vector<double>{.01, .05, .1, .25, .5, .7, .8, .9}))
        cases.emplace_back("bose", z);
    for (double z :
         (quick ? std::vector<double>{1} : std::vector<double>{.01, .05, .1, .25, .5, 1, 2, 5, 10}))
        cases.emplace_back("fermi", z);
    for (const auto &[stat, z] : cases)
        for (int n : {1, 4}) {
            o.statistics = stat;
            o.z = z;
            o.basis = n;
            o.radial = quick ? 12 : (stat == "bose" && z >= .8 ? 40 : 24);
            auto r = solve(o);
            full[key(stat, z, n)] = r;
            records.push_back(r);
            std::cout << "Transport " << stat << " z=" << z << " N=" << n << std::endl;
        }
    for (auto &r : records)
        for (const std::string name : {"viscosity", "conductivity"})
            r[name + "_over_classical_same_order"] =
                r.at(name).get<double>() /
                full.at(key("classical", .1, r.at("basis_size"))).at(name).get<double>();
    write_json(out + "/transport.json", records);
    csv(out + "/transport_same_order.csv", records,
        {"statistics", "fugacity", "basis_size", "viscosity", "conductivity",
         "viscosity_over_classical_same_order", "conductivity_over_classical_same_order"});
    for (const auto &[stat, z] : std::vector<std::pair<std::string, double>>{
             {"classical", .1}, {"bose", .5}, {"bose", .9}, {"fermi", 1}, {"fermi", 10}}) {
        if (quick && (z == .9 || z == 10))
            continue;
        auto r = full.at(key(stat, z, 1));
        auto ref = wu_reference(stat, z, quick ? 32 : 64, quick ? 24 : 40);
        Json b = {{"statistics", stat}, {"fugacity", z}, {"reference", ref}};
        for (const std::string name : {"viscosity", "conductivity"}) {
            double ratio = r.at(name + "_over_classical_first_sonine");
            b[name + "_ratio"] = ratio;
            b[name + "_relative_error"] = ratio / ref.at(name + "_ratio").get<double>() - 1;
        }
        benchmarks.push_back(b);
        std::cout << "Wu benchmark " << stat << " z=" << z << std::endl;
    }
    write_json(out + "/benchmark.json", benchmarks);
    for (const auto &[stat, z] : std::vector<std::pair<std::string, double>>{
             {"classical", .1}, {"bose", .9}, {"fermi", 10}}) {
        o.statistics = stat;
        o.z = z;
        o.basis = 5;
        o.radial = quick ? 16 : 40;
        o.angular = quick ? 12 : 24;
        o.azimuth = quick ? 16 : 28;
        auto fine = solve(o);
        for (int n = 1; n <= 5; ++n) {
            auto tr = truncate(fine, n);
            convergence.push_back({{"kind", "basis"},
                                   {"statistics", stat},
                                   {"fugacity", z},
                                   {"basis_size", n},
                                   {"viscosity", tr.at("viscosity")},
                                   {"conductivity", tr.at("conductivity")}});
        }
        if (!quick)
            for (auto sizes :
                 std::vector<std::array<int, 3>>{{16, 12, 16}, {24, 16, 20}, {40, 24, 28}}) {
                o.basis = 4;
                o.radial = sizes[0];
                o.angular = sizes[1];
                o.azimuth = sizes[2];
                auto r = solve(o);
                convergence.push_back({{"kind", "quadrature"},
                                       {"statistics", stat},
                                       {"fugacity", z},
                                       {"basis_size", 4},
                                       {"radial_order", o.radial},
                                       {"angular_order", o.angular},
                                       {"azimuthal_order", o.azimuth},
                                       {"viscosity", r.at("viscosity")},
                                       {"conductivity", r.at("conductivity")}});
            }
        std::cout << "Convergence " << stat << std::endl;
    }
    write_json(out + "/convergence.json", convergence);
    for (const auto &[stat, z] : std::vector<std::pair<std::string, double>>{
             {"classical", .1}, {"bose", .5}, {"fermi", 1}}) {
        auto r = full.at(key(stat, z, 4));
        double mu0 = r.at("classical_first_sonine_viscosity");
        Json scan = Json::array();
        for (double k : {.05, .1, .2, .4, .8})
            scan.push_back(shear_wave(r, k));
        dynamics.push_back({{"statistics", stat},
                            {"fugacity", z},
                            {"transport", r},
                            {"wave", shear_wave(r, .2, linspace(0, 2.5 / (mu0 * .04), 151))},
                            {"wavelength_scan", scan}});
    }
    write_json(out + "/dynamics.json", dynamics);
    double worst = 0;
    for (const auto &b : benchmarks)
        for (const std::string obs : {"viscosity", "conductivity"})
            worst = std::max(worst, std::abs(b.at(obs + "_relative_error").get<double>()));
    Json manifest = {
        {"implementation", "C++17"},
        {"quick", quick},
        {"benchmark_cases", benchmarks.size()},
        {"max_benchmark_relative_discrepancy", worst},
        {"elapsed_seconds",
         std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()},
        {"scope", "UU Galerkin transport, homogeneous relaxation, projected shear wave"}};
    write_json(out + "/manifest.json", manifest);
    std::ofstream summary(out + "/SUMMARY.md");
    summary << "# C++ reproduction\n\n"
            << (quick ? "Quick development run; not the publication settings."
                      : "Full production quadratures and matched classical/quantum basis orders.")
            << "\n\nMaximum one-profile benchmark discrepancy: " << std::setprecision(8)
            << 100 * worst
            << "%.\n\nThe homogeneous response integral and shear-wave hydrodynamic limit are "
               "internal checks.\n";
    return manifest;
}
Json reproduce_hilker(const std::string &root, bool quick) {
    std::string out = root + "/results/cpp/hilker";
    std::filesystem::create_directories(out);
    auto rows = observations(root + "/references/hilker2022/Figure_4bc_exp.csv"),
         theory = read_table(root + "/references/hilker2022/Figure_4bc_th.csv");
    auto transport = [&](double z, int n, bool fine, bool classical = false) {
        Options o;
        o.statistics = classical ? "classical" : "bose";
        o.z = z;
        o.basis = n;
        o.radial = quick ? 12 : (fine ? 32 : 24);
        o.angular = quick ? 8 : (fine ? 20 : 16);
        o.azimuth = quick ? 10 : (fine ? 24 : 20);
        o.equilibrium = 256;
        return solve(o);
    };
    auto classical = transport(.1, 7, true, true), fixed = transport(.75, 7, true);
    Equilibrium fixed_eq("bose", .75, 256);
    auto gf = bose_functions(.75);
    double cvn = 3.75 * gf[2] / gf[1] - 2.25 * gf[1] / gf[0];
    Json reference = Json::array(), curves = Json::array(), output = Json::array(),
         issues = Json::array();
    double speed_err = 0, damp_err = 0;
    for (const auto &row : theory) {
        double inverse = row.at("L/l_mfp"), qe = std::sqrt(15 * pi) / 12 * .73 * pi / inverse,
               qk = qe * 3 / (4 * .67 * .5);
        CVec roots = sound_modes(qe, qk, .5);
        auto c = roots[0];
        for (int i = 1; i < 3; ++i)
            if (roots[i].real() > c.real())
                c = roots[i];
        Json r = {{"inverse_kn", inverse},
                  {"published_speed", row.at("c_I")},
                  {"published_damping_ratio", row.at("gamma_I/(c_I k)")},
                  {"speed", nullptr},
                  {"damping_ratio", nullptr}};
        if (c.real() > 1e-9) {
            r["speed"] = std::abs(c);
            r["damping_ratio"] = -c.imag() / std::abs(c);
            if (inverse >= 3) {
                speed_err =
                    std::max(speed_err, std::abs(std::abs(c) / row.at("c_I").get<double>() - 1));
                damp_err = std::max(damp_err, std::abs((-c.imag() / std::abs(c)) /
                                                           row.at("gamma_I/(c_I k)").get<double>() -
                                                       1));
            }
        }
        reference.push_back(r);
        for (int n : {1, 7}) {
            auto tr = truncate(fixed, n);
            double eta = tr.at("viscosity").get<double>() /
                         tr.at("classical_first_sonine_viscosity").get<double>() /
                         std::sqrt(fixed_eq.pressure),
                   heat = tr.at("conductivity").get<double>() /
                          tr.at("classical_first_sonine_conductivity").get<double>() /
                          std::sqrt(fixed_eq.pressure);
            auto rr = sound_modes(std::sqrt(15 * pi) / 12 * eta * pi / inverse,
                                  15. / 64 * std::sqrt(15 * pi) * heat / cvn * pi / inverse, .5);
            auto cc = rr[0];
            for (int i = 1; i < 3; ++i)
                if (rr[i].real() > cc.real())
                    cc = rr[i];
            if (cc.real() > 1e-9)
                curves.push_back({{"basis_size", n},
                                  {"inverse_kn", inverse},
                                  {"speed", std::abs(cc)},
                                  {"damping_ratio", -cc.imag() / std::abs(cc)}});
        }
    }
    for (const auto &row : rows) {
        auto state = state_from_observation(row);
        Json item = {
            {"observation", row},
            {"state", state},
            {"hf_speed_relative_to_published",
             state.at("sound_speed").get<double>() / row.at("reference_hf_speed").get<double>() -
                 1}};
        double inconsistency =
            row.at("damping").get<double>() /
                (row.at("speed").get<double>() * pi / row.at("length").get<double>()) -
            row.at("damping_ratio").get<double>();
        if (std::abs(inconsistency) > 1e-12)
            issues.push_back(
                {{"row", row.at("id")},
                 {"issue", "Published damping ratio inconsistent with tabulated length"},
                 {"selected", row.at("hydro_selected")}});
        if (row.at("hydro_selected").get<bool>()) {
            std::cout << "Hilker point " << row.at("id")
                      << ": N=7, two quadratures and input sensitivities" << std::endl;
            double z = state.at("fugacity");
            auto fine = transport(z, 7, true), coarse = transport(z, 7, false),
                 n5 = truncate(fine, 5);
            item["prediction"] = prediction(row, fine);
            item["no_mean_field"] = prediction(row, fine, false);
            item["classical_transport_same_hf"] = prediction(row, classical);
            for (const std::string name : {"viscosity", "conductivity"}) {
                item["convergence"]["basis_5_to_7"][name] =
                    fine.at(name).get<double>() / n5.at(name).get<double>() - 1;
                item["convergence"]["quadrature"][name] =
                    fine.at(name).get<double>() / coarse.at(name).get<double>() - 1;
            }
            item["basis_sequence"] = Json::array();
            for (int n = 1; n <= 7; ++n) {
                auto tr = truncate(fine, n);
                item["basis_sequence"].push_back({{"N", n},
                                                  {"viscosity", tr.at("viscosity")},
                                                  {"conductivity", tr.at("conductivity")}});
            }
            item["sensitivity"] = Json::array();
            for (const std::string name : {"inverse_kn", "scattering_length", "temperature"}) {
                Json changes = Json::array();
                for (int sign : {-1, 1}) {
                    auto shifted = row;
                    shifted[name] =
                        row.at(name).get<double>() + sign * row.at(name + "_error").get<double>();
                    auto st = state_from_observation(shifted);
                    auto r = transport(st.at("fugacity"), 5, false);
                    auto p = prediction(shifted, r);
                    changes.push_back({{"sign", sign},
                                       {"input", shifted.at(name)},
                                       {"fugacity", st.at("fugacity")},
                                       {"speed", p.at("speed")},
                                       {"damping", p.at("damping")}});
                }
                item["sensitivity"].push_back({{"input", name}, {"changes", changes}});
            }
        }
        output.push_back(item);
        write_json(out + "/comparison.json", output);
    }
    // Standing longitudinal mode: density cos(kx), velocity sin(kx), temperature cos(kx).
    Json standing = Json::array();
    for (const auto &item : output)
        if (item.contains("prediction") && item.at("observation").at("id") == 9) {
            auto state = item.at("state"), p = item.at("prediction");
            double l = item.at("observation").at("length"), k = pi / l, c = state.at("sound_speed");
            CMat transform = CMat::Identity(3, 3);
            transform(1, 1) = std::complex<double>(0, -1);
            CMat gen = transform.adjoint() *
                       longitudinal_matrix(state, p.at("viscosity"), p.at("conductivity"), k) *
                       transform * k * c;
            Eigen::ComplexEigenSolver<CMat> eig(gen);
            CVec initial = CVec::Zero(3);
            initial[0] = .01;
            CVec a = eig.eigenvectors().fullPivLu().solve(initial);
            for (double t : linspace(0, .1, 401)) {
                CVec y = eig.eigenvectors() *
                         (a.array() * (t * eig.eigenvalues().array()).exp()).matrix();
                standing.push_back({{"time", t},
                                    {"density_amplitude", y[0].real()},
                                    {"velocity_over_sound_speed", y[1].real()},
                                    {"temperature_amplitude", y[2].real()},
                                    {"length", l}});
            }
        }
    write_json(out + "/standing_wave.json", standing);
    write_json(out + "/reference_curve.json", reference);
    write_json(out + "/fixed_state_galerkin.json", curves);
    Json manifest = {
        {"implementation", "C++17"},
        {"quick", quick},
        {"selection", "L/ell>=3, five observations"},
        {"fitted_parameters", Json::array()},
        {"reference_check_max_relative_errors",
         {{"speed", speed_err}, {"damping_ratio", damp_err}}},
        {"source_consistency_issues", issues},
        {"uncertainty",
         "One-at-a-time sensitivity, not confidence intervals; input covariances unavailable"},
        {"density", "Reconstructed from published L/ell=8*pi*n*a^2*L"}};
    write_json(out + "/manifest.json", manifest);
    return manifest;
}
Json regression_audit(const std::string &root) {
    auto current = read_json(root + "/results/cpp/transport.json");
    std::ifstream f(root + "/tests/reference/transport_same_order.csv");
    require(bool(f), "Missing frozen transport reference");
    std::string line;
    std::getline(f, line);
    std::vector<std::string> header;
    std::stringstream hs(line);
    std::string cell;
    while (std::getline(hs, cell, ','))
        header.push_back(cell);
    Json comparisons = Json::array();
    double max_transport = 0, max_sound = 0;
    while (std::getline(f, line)) {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::map<std::string, std::string> row;
        size_t i = 0;
        while (std::getline(ss, cell, ',')) {
            require(i < header.size(), "Invalid reference CSV");
            row[header[i++]] = cell;
        }
        bool found = false;
        for (const auto &r : current)
            if (r.at("statistics") == row.at("statistics") &&
                r.at("basis_size") == std::stoi(row.at("basis_size")) &&
                std::abs(r.at("fugacity").get<double>() - std::stod(row.at("fugacity"))) < 1e-14) {
                found = true;
                Json diff = {{"statistics", r.at("statistics")},
                             {"fugacity", r.at("fugacity")},
                             {"basis_size", r.at("basis_size")}};
                for (const std::string name : {"viscosity", "conductivity"}) {
                    double error = std::abs(r.at(name).get<double>() / std::stod(row.at(name)) - 1);
                    diff[name + "_relative_difference"] = error;
                    max_transport = std::max(max_transport, error);
                }
                comparisons.push_back(diff);
                break;
            }
        require(found, "Full C++ reproduction is required before audit");
    }
    auto old = read_json(root + "/tests/reference/hilker_comparison.json"),
         fresh = read_json(root + "/results/cpp/hilker/comparison.json");
    int sound_count = 0;
    for (const auto &r : old)
        if (r.contains("prediction")) {
            bool found = false;
            for (const auto &s : fresh)
                if (s.at("observation").at("id") == r.at("observation").at("id")) {
                    for (const std::string name : {"speed", "damping", "viscosity", "conductivity"})
                        max_sound = std::max(
                            max_sound, std::abs(s.at("prediction").at(name).get<double>() /
                                                    r.at("prediction").at(name).get<double>() -
                                                1));
                    found = true;
                    ++sound_count;
                }
            require(found, "Missing C++ sound prediction");
        }
    Json report = {
        {"transport_cases", comparisons.size()},
        {"max_transport_relative_difference", max_transport},
        {"hilker_cases", sound_count},
        {"max_hilker_prediction_relative_difference", max_sound},
        {"transport", comparisons},
        {"tolerance", 2e-8},
        {"passed", max_transport < 2e-8 && max_sound < 2e-8},
        {"interpretation",
         "Regression consistency with frozen prototype results; not independent physical validation"}};
    write_json(root + "/results/cpp/regression_audit.json", report);
    require(report.at("passed").get<bool>(), "Regression tolerance exceeded");
    return report;
}
} // namespace qt
