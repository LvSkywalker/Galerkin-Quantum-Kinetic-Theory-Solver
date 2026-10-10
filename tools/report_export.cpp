// Export already computed C++ results to PGFPlots data and LaTeX tables.
// No collision calculation is repeated and no Python runtime is used.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
using Json = nlohmann::json;
Json read(const std::string &p) {
    std::ifstream f(p);
    if (!f)
        throw std::runtime_error("Missing " + p);
    Json j;
    f >> j;
    return j;
}
std::ofstream output(const std::string &p) {
    std::ofstream f(p);
    if (!f)
        throw std::runtime_error("Cannot write " + p);
    f << std::setprecision(16);
    return f;
}
int main(int argc, char **argv) {
    try {
        if (argc != 3)
            throw std::runtime_error(
                "Usage: qtransport_report_export PROJECT_ROOT OUTPUT_DIRECTORY");
        std::string root = argv[1], dest = argv[2], data = dest + "/data/";
        std::filesystem::create_directories(data);
        auto transport = read(root + "/results/cpp/transport.json"),
             conv = read(root + "/results/cpp/convergence.json"),
             bench = read(root + "/results/cpp/benchmark.json"),
             dyn = read(root + "/results/cpp/dynamics.json"),
             hilker = read(root + "/results/cpp/hilker/comparison.json");
        for (std::string stat : {"bose", "fermi"})
            for (int n : {1, 4}) {
                auto f = output(data + stat + "_N" + std::to_string(n) + ".dat");
                f << "z mu lambda\n";
                for (const auto &r : transport)
                    if (r["statistics"] == stat && r["basis_size"] == n)
                        f << r["fugacity"] << ' ' << r["viscosity_over_classical_same_order"] << ' '
                          << r["conductivity_over_classical_same_order"] << '\n';
            }
        auto ct = output(data + "convergence_table.tex");
        ct << std::fixed << std::setprecision(5);
        for (std::string stat : {"classical", "bose", "fermi"}) {
            Json n4, n5, coarse, fine;
            for (const auto &r : conv)
                if (r["statistics"] == stat) {
                    if (r["kind"] == "basis" && r["basis_size"] == 4)
                        n4 = r;
                    if (r["kind"] == "basis" && r["basis_size"] == 5)
                        n5 = r;
                    if (r["kind"] == "quadrature" && r["radial_order"] == 24)
                        coarse = r;
                    if (r["kind"] == "quadrature" && r["radial_order"] == 40)
                        fine = r;
                }
            auto f = output(data + "convergence_" + stat + ".dat");
            f << "N mu lambda\n";
            for (const auto &r : conv)
                if (r["kind"] == "basis" && r["statistics"] == stat)
                    f << r["basis_size"] << ' '
                      << 100 * (r["viscosity"].get<double>() / n5["viscosity"].get<double>() - 1)
                      << ' '
                      << 100 * (r["conductivity"].get<double>() / n5["conductivity"].get<double>() -
                                1)
                      << '\n';
            ct << stat << " & " << (stat == "classical" ? "--" : n4["fugacity"].dump());
            for (const std::string obs : {"viscosity", "conductivity"})
                ct << " & " << 100 * std::abs(n5[obs].get<double>() / n4[obs].get<double>() - 1);
            for (const std::string obs : {"viscosity", "conductivity"})
                ct << " & "
                   << 100 * std::abs(fine[obs].get<double>() / coarse[obs].get<double>() - 1);
            ct << " \\\\\n";
        }
        auto bt = output(data + "benchmark_table.tex");
        bt << std::scientific << std::setprecision(2);
        for (const auto &r : bench)
            bt << r["statistics"].get<std::string>() << " & "
               << (r["statistics"] == "classical" ? "--" : r["fugacity"].dump()) << " & $"
               << std::abs(r["viscosity_relative_error"].get<double>()) << "$ & $"
               << std::abs(r["conductivity_relative_error"].get<double>()) << "$ \\\\\n";
        auto tr = output(data + "transport_table.tex");
        tr << std::fixed << std::setprecision(6);
        for (const auto &r : transport)
            if (r["basis_size"] == 4 &&
                ((r["statistics"] == "classical") ||
                 (r["statistics"] == "bose" && (r["fugacity"] == .5 || r["fugacity"] == .9)) ||
                 (r["statistics"] == "fermi" && (r["fugacity"] == 1. || r["fugacity"] == 10.))))
                tr << r["statistics"].get<std::string>() << " & "
                   << (r["statistics"] == "classical" ? "--" : r["fugacity"].dump()) << " & "
                   << r["viscosity"].get<double>() << " & " << r["conductivity"].get<double>()
                   << " & " << r["viscosity_over_classical_same_order"].get<double>() << " & "
                   << r["conductivity_over_classical_same_order"].get<double>() << " \\\\\n";
        auto rates = output(data + "relaxation_table.tex");
        rates << std::fixed << std::setprecision(6);
        auto waves = output(data + "wave_table.tex");
        waves << std::fixed << std::setprecision(6);
        for (const auto &r : dyn) {
            std::string stat = r["statistics"];
            const auto &t = r["transport"];
            const auto &h = t["heat_relaxation"];
            const auto &s = t["shear_relaxation"];
            auto f = output(data + "relaxation_" + stat + ".dat");
            f << "t heat shear Hheat Hshear\n";
            for (size_t i = 0; i < h["time"].size(); ++i)
                f << h["time"][i] << ' ' << h["normalized_current"][i] << ' '
                  << s["normalized_current"][i] << ' ' << h["normalized_quadratic_entropy"][i]
                  << ' ' << s["normalized_quadratic_entropy"][i] << '\n';
            rates << stat << " & " << (stat == "classical" ? "--" : r["fugacity"].dump()) << " & "
                  << h["rates"][0].get<double>() << " & " << s["rates"][0].get<double>()
                  << " \\\\\n";
            auto g = output(data + "wave_" + stat + ".dat");
            g << "t amplitude hydro\n";
            auto w = r["wave"];
            double scale = t["classical_first_sonine_viscosity"].get<double>() * .04;
            for (size_t i = 0; i < w["time"].size(); ++i)
                g << scale * w["time"][i].get<double>() << ' ' << w["velocity_amplitude"][i] << ' '
                  << w["navier_stokes_amplitude"][i] << '\n';
            auto scan = output(data + "wave_scan_" + stat + ".dat");
            scan << "k error\n";
            for (const auto &entry : r["wavelength_scan"]) {
                scan << entry["wave_number"] << ' '
                     << 100 * entry["relative_difference_from_transport_viscosity"].get<double>()
                     << '\n';
                if (entry["wave_number"] == .05)
                    waves << stat << " & " << t["viscosity"].get<double>() << " & "
                          << entry["viscosity_from_slow_mode"].get<double>() << " & "
                          << 100 *
                                 entry["relative_difference_from_transport_viscosity"].get<double>()
                          << " \\\\\n";
            }
        }
        auto hf = output(data + "hilker.dat"), ht = output(data + "hilker_table.tex");
        hf << "index inverse_kn length cdata cerr cmodel csens cclassical gdata gerr gmodel gsens gclassical\n";
        ht << std::fixed;
        int index = 0;
        constexpr double pi = 3.14159265358979323846;
        for (const auto &item : hilker)
            if (item.contains("prediction")) {
                auto r = item["observation"], p = item["prediction"];
                double cs = 0, gs = 0;
                for (const auto &sensitivity : item["sensitivity"])
                    for (const auto &change : sensitivity["changes"]) {
                        cs = std::max(
                            cs, std::abs(change["speed"].get<double>() - p["speed"].get<double>()));
                        gs = std::max(gs, std::abs(change["damping"].get<double>() -
                                                   p["damping"].get<double>()));
                    }
                hf << ++index << ' ' << r["inverse_kn"].get<double>() << ' '
                   << std::lround(1e6 * r["length"].get<double>()) << ' '
                   << 1e3 * r["speed"].get<double>() << ' '
                   << 1e3 * r["speed_error"].get<double>() << ' ' << 1e3 * p["speed"].get<double>()
                   << ' ' << 1e3 * cs << ' '
                   << 1e3 * item.at("classical_transport_same_hf").at("speed").get<double>() << ' '
                   << r["damping"].get<double>() / (2 * pi) << ' '
                   << r["damping_error"].get<double>() / (2 * pi) << ' '
                   << p["damping"].get<double>() / (2 * pi) << ' ' << gs / (2 * pi) << ' '
                   << item.at("classical_transport_same_hf").at("damping").get<double>() / (2 * pi)
                   << '\n';
                ht << index << " & " << std::setprecision(0) << 1e6 * r["length"].get<double>()
                   << " & " << r["raw"]["a (aB)"].get<double>() << " & " << std::setprecision(3)
                   << 1e3 * r["speed"].get<double>() << " & " << 1e3 * p["speed"].get<double>()
                   << " & " << std::setprecision(2) << r["damping"].get<double>() / (2 * pi)
                   << " & " << p["damping"].get<double>() / (2 * pi) << " \\\\\n";
            }
        Json metrics = {{"selected_observations", index},
                        {"selection", "L/ell >= 3"},
                        {"normalization", "Each residual divided by its reported measurement uncertainty"},
                        {"control", "Classical transport with the SAME Hartree-Fock thermodynamics"},
                        {"interpretation", "Descriptive residuals; input uncertainty and covariance not propagated"}};
        auto mt = output(data + "hilker_metrics_table.tex");
        mt << std::fixed << std::setprecision(3);
        for (const std::string model : {"prediction", "classical_transport_same_hf"}) {
            const std::string name = model == "prediction" ? "quantum_transport" : "classical_transport";
            for (const std::string obs : {"speed", "damping"}) {
                double sumsq = 0, minrel = std::numeric_limits<double>::infinity(), maxrel = -minrel;
                size_t count = 0;
                Json residuals = Json::array();
                for (const auto &item : hilker) {
                    if (!item.contains("prediction"))
                        continue;
                    const auto &measurement = item.at("observation");
                    double value = measurement.at(obs), error = measurement.at(obs + "_error"),
                           predicted = item.at(model).at(obs);
                    if (!std::isfinite(value) || value <= 0 || !std::isfinite(error) || error <= 0 ||
                        !std::isfinite(predicted))
                        throw std::runtime_error("Invalid sound comparison value or uncertainty");
                    double residual = (predicted - value) / error, relative = predicted / value - 1;
                    sumsq += residual * residual;
                    minrel = std::min(minrel, relative);
                    maxrel = std::max(maxrel, relative);
                    residuals.push_back({{"observation_id", measurement.at("id")},
                                         {"standardized_residual", residual}, {"relative_residual", relative}});
                    ++count;
                }
                if (!count)
                    throw std::runtime_error("No selected sound observations");
                metrics[name][obs] = {{"rms_standardized_residual", std::sqrt(sumsq / count)},
                                      {"min_relative_residual", minrel}, {"max_relative_residual", maxrel},
                                      {"max_absolute_relative_residual", std::max(std::abs(minrel), std::abs(maxrel))},
                                      {"residuals", residuals}};
            }
            mt << (model == "prediction" ? "Quantum transport" : "Classical transport")
               << " & " << metrics[name]["speed"]["rms_standardized_residual"].get<double>()
               << " & " << metrics[name]["damping"]["rms_standardized_residual"].get<double>()
               << " \\\\\n";
        }
        auto mj = output(data + "hilker_metrics.json");
        mj << metrics.dump(2) << '\n';
        std::cout << "Exported C++ results and report tables to " << dest << '\n';
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
