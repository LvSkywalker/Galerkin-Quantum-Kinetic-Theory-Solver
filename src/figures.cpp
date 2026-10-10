#include "qtransport/core.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace qt {
namespace {
struct Curve {
    std::string name, color;
    std::vector<double> x, y, error;
};
std::string number(double x) {
    std::ostringstream s;
    s << std::setprecision(3) << x;
    return s.str();
}
struct Svg {
    std::ofstream f;
    int width, height;
    Svg(const std::string &path, int w, int h, const std::string &title)
        : f(path), width(w), height(h) {
        if (!f)
            throw std::runtime_error("Cannot write SVG");
        f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
          << "\" viewBox=\"0 0 " << w << ' ' << h
          << "\"><rect width=\"100%\" height=\"100%\" "
             "fill=\"white\"/"
             "><style>text{font-family:Arial,sans-serif;fill:#243044;font-size:13px}.title{font-"
             "size:21px;font-weight:bold}.panel{font-size:16px;font-weight:bold}</style>";
        text(w / 2, 30, title, "title", "middle");
    }
    ~Svg() {
        f << "</svg>\n";
    }
    void text(double x, double y, const std::string &s, const std::string &cls = "",
              const std::string &anchor = "start") {
        f << "<text x=\"" << x << "\" y=\"" << y << "\" class=\"" << cls << "\" text-anchor=\""
          << anchor << "\">" << s << "</text>";
    }
    void line(double x, double y, double xx, double yy, const std::string &color, double w = 1) {
        f << "<path d=\"M" << x << ',' << y << " L" << xx << ',' << yy << "\" stroke=\"" << color
          << "\" stroke-width=\"" << w << "\"/>";
    }
    void panel(double left, double top, double w, double h, const std::string &title,
               const std::string &xlabel, const std::string &ylabel,
               const std::vector<Curve> &curves, bool logy = false) {
        double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
        for (const auto &c : curves)
            for (size_t i = 0; i < c.x.size(); ++i) {
                if (!std::isfinite(c.y[i]) || (logy && c.y[i] <= 0))
                    continue;
                xmin = std::min(xmin, c.x[i]);
                xmax = std::max(xmax, c.x[i]);
                double error = c.error.empty() ? 0 : c.error[i];
                ymin = std::min(ymin, logy ? std::log10(c.y[i]) : c.y[i] - error);
                ymax = std::max(ymax, logy ? std::log10(c.y[i]) : c.y[i] + error);
            }
        if (!(xmax > xmin)) {
            xmin = 0;
            xmax = 1;
        }
        if (!(ymax > ymin)) {
            ymin -= .5;
            ymax += .5;
        }
        double pad = .04 * (ymax - ymin);
        ymin -= pad;
        ymax += pad;
        if (logy) {
            ymin = std::floor(ymin);
            ymax = std::ceil(ymax);
        }
        auto px = [&](double x) { return left + w * (x - xmin) / (xmax - xmin); };
        auto py = [&](double y) { return top + h - h * (y - ymin) / (ymax - ymin); };
        text(left, top - 17, title, "panel");
        text(left, top - 1, ylabel);
        for (int j = 0; j <= 4; ++j) {
            double x = xmin + (xmax - xmin) * j / 4, y = ymin + (ymax - ymin) * j / 4;
            line(px(x), top, px(x), top + h, "#e3e7ed");
            text(px(x), top + h + 20, number(x), "", "middle");
            if (!logy) {
                line(left, py(y), left + w, py(y), "#e3e7ed");
                text(left - 7, py(y) + 4, number(y), "", "end");
            }
        }
        if (logy) {
            int step = std::max(1, int(std::ceil((ymax - ymin) / 5)));
            for (int y = int(ymax); y >= int(ymin); y -= step) {
                line(left, py(y), left + w, py(y), "#e3e7ed");
                text(left - 7, py(y) + 4, "10^" + std::to_string(y), "", "end");
            }
        }
        line(left, top, left, top + h, "#243044");
        line(left, top + h, left + w, top + h, "#243044");
        text(left + w / 2, top + h + 42, xlabel, "", "middle");
        int legend = 0;
        for (const auto &c : curves) {
            f << "<polyline fill=\"none\" stroke=\"" << c.color
              << "\" stroke-width=\"2\" points=\"";
            for (size_t i = 0; i < c.x.size(); ++i)
                if (std::isfinite(c.y[i]) && (!logy || c.y[i] > 0))
                    f << px(c.x[i]) << ',' << py(logy ? std::log10(c.y[i]) : c.y[i]) << ' ';
            f << "\"/>";
            if (!c.error.empty())
                for (size_t i = 0; i < c.x.size(); ++i) {
                    double x = px(c.x[i]), y = py(c.y[i]);
                    line(x, py(c.y[i] - c.error[i]), x, py(c.y[i] + c.error[i]), c.color);
                    f << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"3\" fill=\"" << c.color
                      << "\"/>";
                }
            double ly = top + h - 12 - 18 * legend++;
            line(left + 12, ly - 4, left + 34, ly - 4, c.color, 2);
            text(left + 40, ly, c.name);
        }
    }
};
std::string color(const std::string &stat) {
    return stat == "bose" ? "#087f8c" : stat == "fermi" ? "#b24c63" : "#626b78";
}
std::vector<double> values(const Json &j) {
    return j.get<std::vector<double>>();
}
} // namespace
void make_figures(const std::string &root) {
    std::string out = root + "/results/cpp";
    auto dynamics = read_json(out + "/dynamics.json"), records = read_json(out + "/transport.json");
    {
        Svg svg(out + "/homogeneous_relaxation.svg", 1100, 800,
                "Homogeneous linearized relaxation — C++ Galerkin solver");
        svg.text(550, 54, "Same kernel B = d² |v-v₁| / 4; N = 4 radial profiles; c(0) = H⁻¹d", "",
                 "middle");
        for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 2; ++col) {
                std::vector<Curve> curves;
                std::string sector = col ? "shear" : "heat";
                for (const auto &e : dynamics) {
                    auto r = e.at("transport");
                    std::string stat = r.at("statistics"), label = stat;
                    if (stat != "classical")
                        label += " z=" + number(r.at("fugacity"));
                    auto data = r.at(sector + "_relaxation");
                    curves.push_back({label,
                                      color(stat),
                                      values(data.at("time")),
                                      values(data.at(row ? "normalized_quadratic_entropy"
                                                         : "normalized_current")),
                                      {}});
                }
                svg.panel(90 + 535 * col, 110 + 350 * row, 430, 250,
                          col ? "Shear-stress perturbation" : "Heat-flux perturbation",
                          "t / t₀; t₀ = (n d² √(kBT/m))⁻¹",
                          row ? "Quadratic relative entropy / initial value"
                              : "Flux / initial value",
                          curves, true);
            }
    }
    {
        Svg svg(out + "/transport.svg", 1100, 780,
                "Quantum transport — classical reference at the same basis order");
        for (int col = 0; col < 2; ++col)
            for (int row = 0; row < 2; ++row) {
                std::string stat = col ? "fermi" : "bose", obs = row ? "conductivity" : "viscosity";
                std::vector<Curve> curves;
                for (int n : {1, 4}) {
                    Curve c{"N=" + std::to_string(n), n == 1 ? "#8795a9" : color(stat), {}, {}, {}};
                    for (const auto &r : records)
                        if (r.at("statistics") == stat && r.at("basis_size") == n) {
                            c.x.push_back(r.at("fugacity"));
                            c.y.push_back(r.at(obs + "_over_classical_same_order"));
                        }
                    curves.push_back(c);
                }
                svg.panel(80 + 540 * col, 100 + 340 * row, 440, 240, stat + ": " + obs,
                          "Fugacity z", "Quantum / classical coefficient", curves);
            }
    }
    {
        Svg svg(out + "/shear_wave.svg", 1100, 430,
                "Kinetic shear-wave decay and hydrodynamic consistency");
        std::vector<Curve> decay, scan;
        for (const auto &e : dynamics) {
            std::string stat = e.at("statistics");
            auto wave = e.at("wave");
            Curve c{stat,
                    color(stat),
                    values(wave.at("time")),
                    values(wave.at("velocity_amplitude")),
                    {}};
            double scale =
                e.at("transport").at("classical_first_sonine_viscosity").get<double>() * .04;
            for (auto &x : c.x)
                x *= scale;
            decay.push_back(c);
            Curve b{stat, color(stat), {}, {}, {}};
            for (const auto &s : e.at("wavelength_scan")) {
                b.x.push_back(s.at("wave_number"));
                b.y.push_back(100 *
                              s.at("relative_difference_from_transport_viscosity").get<double>());
            }
            scan.push_back(b);
        }
        svg.panel(80, 100, 440, 240, "Velocity amplitude", "μclassical,1 k² t", "u(t) / u(0)",
                  decay);
        svg.panel(620, 100, 440, 240, "Viscosity inferred from slow mode", "Wave number k",
                  "Difference from steady viscosity (%)", scan);
    }
    {
        Svg svg(out + "/shear_wave_animation.svg", 850, 420,
                "Projected kinetic shear wave — amplitude evolves from collisions");
        svg.text(425, 55, "Time runs over the saved interval and repeats; x is spatial phase", "",
                 "middle");
        for (const auto &e : dynamics) {
            std::string stat = e.at("statistics");
            auto wave = e.at("wave");
            svg.f << "<polyline fill=\"none\" stroke=\"" << color(stat)
                  << "\" stroke-width=\"2\"><animate attributeName=\"points\" dur=\"12s\" "
                     "repeatCount=\"indefinite\" values=\"";
            int frame = 0;
            for (const auto &a : wave.at("velocity_amplitude")) {
                if (frame++)
                    svg.f << ';';
                for (int i = 0; i <= 100; ++i)
                    svg.f << 60 + 7.3 * i << ','
                          << 215 - 120 * a.get<double>() * std::cos(2 * pi * i / 100) << ' ';
            }
            svg.f << "\"/></polyline>";
        }
        svg.line(60, 215, 790, 215, "#adb5c2");
        svg.text(425, 390, "0 ≤ ky ≤ 2π; curves: classical (grey), Bose (teal), Fermi (red)", "",
                 "middle");
    }
    if (std::filesystem::exists(out + "/convergence.json")) {
        auto data = read_json(out + "/convergence.json");
        Svg svg(out + "/basis_convergence.svg", 1100, 430,
                "Radial-basis refinement at fixed quadrature");
        for (int col = 0; col < 2; ++col) {
            std::string obs = col ? "conductivity" : "viscosity";
            std::vector<Curve> curves;
            for (std::string stat : {"classical", "bose", "fermi"}) {
                Curve c{stat, color(stat), {}, {}, {}};
                double final = 0;
                for (const auto &r : data)
                    if (r.at("kind") == "basis" && r.at("statistics") == stat &&
                        r.at("basis_size") == 5)
                        final = r.at(obs);
                for (const auto &r : data)
                    if (r.at("kind") == "basis" && r.at("statistics") == stat) {
                        c.x.push_back(r.at("basis_size"));
                        c.y.push_back(100 * (r.at(obs).get<double>() / final - 1));
                    }
                curves.push_back(c);
            }
            svg.panel(80 + 540 * col, 100, 440, 240, obs, "Number of radial profiles",
                      "Difference from N=5 (%)", curves);
        }
    }
    if (std::filesystem::exists(out + "/benchmark.json")) {
        auto data = read_json(out + "/benchmark.json");
        Svg svg(out + "/benchmark.svg", 1100, 430,
                "Independent one-profile benchmark — Wu Eq. (50)");
        for (int col = 0; col < 2; ++col) {
            std::string name = col ? "conductivity" : "viscosity";
            Curve c{name, col ? "#b24c63" : "#087f8c", {}, {}, {}};
            for (const auto &r : data) {
                c.x.push_back(c.x.size() + 1);
                c.y.push_back(
                    std::max(1e-16, std::abs(r.at(name + "_relative_error").get<double>())));
            }
            svg.panel(80 + 540 * col, 100, 440, 240, name, "Benchmark case index",
                      "Absolute relative discrepancy", {c}, true);
        }
    }
    if (std::filesystem::exists(out + "/fourier.json")) {
        auto data = read_json(out + "/fourier.json");
        Svg svg(out + "/fourier.svg", 1100, 430,
                "Velocity Fourier representation — off-grid error, not a collision solver");
        for (int col = 0; col < 2; ++col) {
            std::string sector = col ? "shear" : "equilibrium";
            std::vector<Curve> curves;
            for (std::string stat : {"bose", "fermi"})
                for (double width : {6., 8.}) {
                    Curve c{stat + " L=" + number(width),
                            width == 6 ? color(stat) : (stat == "bose" ? "#51b7ab" : "#d69cb0"),
                            {},
                            {},
                            {}};
                    for (const auto &r : data)
                        if (r.at("statistics") == stat && r.at("half_width") == width) {
                            c.x.push_back(r.at("points_per_axis"));
                            c.y.push_back(r.at(sector + "_offgrid_relative_l2"));
                        }
                    curves.push_back(c);
                }
            svg.panel(80 + 540 * col, 100, 440, 240, sector, "Grid points per velocity axis",
                      "Relative L2 interpolation error", curves, true);
        }
    }
    if (std::filesystem::exists(out + "/hilker/comparison.json")) {
        auto data = read_json(out + "/hilker/comparison.json");
        Svg svg(out + "/hilker/comparison.svg", 1100, 450,
                "Normal Bose sound — five selected Hilker observations, no fitted parameters");
        for (int col = 0; col < 2; ++col) {
            std::vector<Curve> curves{{"Experiment", "#243044", {}, {}, {}},
                                      {"C++ Galerkin + HF + NSF", "#c97916", {}, {}, {}}};
            for (const auto &item : data)
                if (item.contains("prediction")) {
                    const auto &r = item.at("observation");
                    const auto &p = item.at("prediction");
                    std::string obs = col ? "damping" : "speed";
                    double scale = col ? 1 / (2 * pi) : 1e3, x = curves[0].x.size() + 1;
                    for (auto &c : curves)
                        c.x.push_back(x);
                    curves[0].y.push_back(scale * r.at(obs).get<double>());
                    curves[0].error.push_back(scale * r.at(obs + "_error").get<double>());
                    curves[1].y.push_back(scale * p.at(obs).get<double>());
                    double sensitivity = 0;
                    for (const auto &s : item.at("sensitivity"))
                        for (const auto &change : s.at("changes"))
                            sensitivity = std::max(sensitivity,
                                                   std::abs(scale * (change.at(obs).get<double>() -
                                                                     p.at(obs).get<double>())));
                    curves[1].error.push_back(sensitivity);
                }
            svg.panel(80 + 540 * col, 105, 440, 250, col ? "Sound damping" : "Sound speed",
                      "Selected observation index", col ? "γ / (2π) [Hz]" : "|ω| / k [mm/s]",
                      curves);
        }
        svg.text(550, 430,
                 "Orange ranges: one-input-at-a-time sensitivity, not confidence intervals.", "",
                 "middle");
    }
    if (std::filesystem::exists(out + "/hilker/standing_wave.json")) {
        auto data = read_json(out + "/hilker/standing_wave.json");
        Svg svg(out + "/hilker/standing_wave.svg", 1100, 500,
                "Longitudinal standing sound wave — linearized Navier–Stokes–Fourier");
        std::vector<Curve> curves;
        int j = 0;
        for (std::string key :
             {"density_amplitude", "velocity_over_sound_speed", "temperature_amplitude"}) {
            Curve c{
                key, std::vector<std::string>{"#087f8c", "#b24c63", "#626b78"}[j++], {}, {}, {}};
            for (const auto &r : data) {
                c.x.push_back(1000 * r.at("time").get<double>());
                c.y.push_back(r.at(key));
            }
            curves.push_back(c);
        }
        svg.panel(95, 100, 950, 300, "Mode amplitudes; density initially perturbed by 1%",
                  "Time [ms]", "Dimensionless field amplitude", curves);
    }
    std::ofstream html(out + "/index.html");
    html << "<!doctype html><meta charset='utf-8'><title>Quantum Gas Transport — C++ "
            "results</title><style>body{font-family:system-ui;max-width:1100px;margin:30px "
            "auto;padding:20px;background:#f4f6f8}img,object{width:100%;background:white;margin-"
            "bottom:24px}p{line-height:1.6}</style><h1>Quantum Gas Transport — C++ "
            "results</h1><p>Deterministic Galerkin collision integrals, conservation constraints "
            "and linearized dynamics. Build settings and numerical checks are documented in the "
            "repository technical note.</p>";
    for (std::string name :
         {"transport", "basis_convergence", "benchmark", "homogeneous_relaxation", "shear_wave",
          "fourier", "shear_wave_animation", "hilker/comparison", "hilker/standing_wave"})
        if (std::filesystem::exists(out + "/" + name + ".svg"))
            html << "<object type='image/svg+xml' data='" << name << ".svg'></object>";
}
} // namespace qt
