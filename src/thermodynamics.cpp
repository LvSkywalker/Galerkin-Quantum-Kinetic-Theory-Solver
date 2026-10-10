#include "qtransport/core.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace qt {
std::array<double, 3> bose_functions(double z) {
    require(std::isfinite(z) && z > 0 && z < 1, "Normal Bose series requires 0<z<1");
    std::array<double, 3> sums{0, 0, 0};
    double power = z;
    for (int j = 1; j <= 500000; ++j) {
        double term = power / std::sqrt(double(j));
        sums[0] += term;
        sums[1] += term / j;
        sums[2] += term / (double(j) * j);
        if (term * z / (1 - z) < 2e-14 * sums[2])
            return sums;
        power *= z;
    }
    throw std::runtime_error("Bose series too close to criticality");
}
Json normal_state(double n, double t, double a, double mass) {
    require(std::isfinite(n) && std::isfinite(t) && std::isfinite(a) && std::isfinite(mass) &&
                n > 0 && t > 0 && a >= 0 && mass > 0,
            "Invalid thermodynamic state");
    double wavelength = std::sqrt(2 * pi * hbar * hbar / (mass * kb * t)),
           phase = n * std::pow(wavelength, 3);
    require(phase > 0 && phase < 2.612375348685488, "State outside normal Bose branch");
    double lo = 0, hi = .9;
    while (bose_functions(hi)[1] < phase)
        hi = 1 - (1 - hi) / 2;
    for (int i = 0; i < 48; ++i) {
        double mid = (lo + hi) / 2;
        if (bose_functions(mid)[1] < phase)
            lo = mid;
        else
            hi = mid;
    }
    double z = (lo + hi) / 2;
    auto g = bose_functions(z);
    double interaction = 4 * pi * hbar * hbar * a / mass;
    double pkin = n * kb * t * g[2] / g[1], pn = kb * t * g[1] / g[0] + 2 * interaction * n,
           pt = n * kb * (2.5 * g[2] / g[1] - 1.5 * g[1] / g[0]), cv = 1.5 * pt,
           cp = cv + t * pt * pt / (n * pn), cad = std::sqrt(cp / cv * pn / mass);
    return {{"n", n},
            {"temperature", t},
            {"scattering_length", a},
            {"mass", mass},
            {"rho", mass * n},
            {"fugacity", z},
            {"thermal_wavelength", wavelength},
            {"phase_density", g[1]},
            {"kinetic_pressure", pkin},
            {"pressure", pkin + interaction * n * n},
            {"energy_density", 1.5 * pkin + interaction * n * n},
            {"chemical_potential", kb * t * std::log(z) + 2 * interaction * n},
            {"pressure_density_derivative", pn},
            {"pressure_temperature_derivative", pt},
            {"cv", cv},
            {"cp", cp},
            {"sound_speed", cad},
            {"mean_field_parameter", interaction * n / (kb * t)},
            {"gas_parameter", n * a * a * a},
            {"thermal_ka", std::sqrt(mass * kb * t) * a / hbar}};
}
std::pair<double, double> dimensional_transport(const Json &r, double t, double a, double mass) {
    require(r.value("kernel_gamma", 1.) == 1 && a > 0 && t > 0 && mass > 0,
            "SI mapping requires gamma=1, positive a,T,m");
    double c = r.value("kernel_strength", .25);
    return {r.at("viscosity").get<double>() * c * std::sqrt(mass * kb * t) / (2 * a * a),
            r.at("conductivity").get<double>() * c * kb * std::sqrt(kb * t / mass) / (2 * a * a)};
}
namespace {
std::vector<std::string> split_csv(const std::string &line) {
    std::vector<std::string> cells;
    std::string cell;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                cell += '"';
                ++i;
            } else
                quoted = !quoted;
        } else if (c == ',' && !quoted) {
            cells.push_back(cell);
            cell.clear();
        } else if (c != '\r')
            cell += c;
    }
    cells.push_back(cell);
    return cells;
}
} // namespace
Json read_table(const std::string &path) {
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Missing Hilker CSV: " + path);
    std::string line;
    std::vector<std::string> header;
    Json rows = Json::array();
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line == "\r")
            continue;
        if (header.empty()) {
            header = split_csv(line);
            continue;
        }
        auto cells = split_csv(line);
        require(cells.size() == header.size(), "Invalid experimental CSV row");
        Json raw;
        for (size_t i = 0; i < header.size(); ++i)
            raw[header[i]] = std::stod(cells[i]);
        rows.push_back(raw);
    }
    return rows;
}
Json observations(const std::string &path) {
    Json rows = Json::array();
    for (const auto &raw : read_table(path)) {
        auto v = [&](const std::string &key) { return raw.at(key).get<double>(); };
        rows.push_back({{"id", rows.size()},
                        {"raw", raw},
                        {"length", v("box L (um)") * 1e-6},
                        {"scattering_length", v("a (aB)") * a0},
                        {"scattering_length_error", v("a uncert. (aB)") * a0},
                        {"inverse_kn", v("L/l_mfp")},
                        {"inverse_kn_error", v("L/l_mfp uncert.")},
                        {"temperature", v("temperature (nK)") * 1e-9},
                        {"temperature_error", v("temperature uncert. (nK)") * 1e-9},
                        {"speed", v("speed c_I (mm/s)") * 1e-3},
                        {"speed_error", v("c_I uncert. (mm/s)") * 1e-3},
                        {"reference_hf_speed", v("c_I^HF (mm/s)") * 1e-3},
                        {"speed_ratio", v("c_I/c_I^HF")},
                        {"speed_ratio_error", v("c_I/c_I^HF uncert.")},
                        {"damping", v("damping gamma_I (2pi Hz)") * 2 * pi},
                        {"damping_error", v("gamma_I uncert.") * 2 * pi},
                        {"damping_ratio", v("gamma_I/(c_I k)")},
                        {"damping_ratio_error", v("gamma_I/(c_I k) uncert.")},
                        {"hydro_selected", v("L/l_mfp") >= 3}});
    }
    return rows;
}
Json state_from_observation(const Json &r, bool mean_field) {
    double a = r.at("scattering_length"), l = r.at("length"), kn = r.at("inverse_kn");
    return normal_state(kn / (8 * pi * a * a * l), r.at("temperature"), mean_field ? a : 0);
}
} // namespace qt
