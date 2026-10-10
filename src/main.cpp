#include "qtransport/core.hpp"
#include <iostream>
#include <map>
#include <stdexcept>

int main(int argc, char **argv) {
    try {
        if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
            std::cout
                << "Quantum Gas Transport -- C++17\n\n"
                   "Usage: qtransport COMMAND [options]\n\n"
                   "Commands:\n"
                   "  transport   Galerkin coefficients and homogeneous linear relaxation\n"
                   "  reference   Wu Eq. (50) one-profile reference integrals\n"
                   "  reproduce   Production scan, benchmark and convergence studies\n"
                   "  figures     SVG figures and gallery from saved results\n"
                   "  test        Automated physics and regression checks\n"
                   "  audit       Compare production outputs with frozen prototype data\n"
                   "  wave        Projected transverse kinetic mode\n"
                   "  hilker      Optional normal-Bose sound comparison\n\n"
                   "Options: --statistics classical|bose|fermi --fugacity Z --basis N\n"
                   " --radial N --angular N --azimuth N --equilibrium N --gamma G --strength C\n"
                   " --output FILE --input FILE --root DIRECTORY --quick --k K\n"
                   "Defaults: classical, Z=0.1, basis=2, radial=16, angular=12, azimuth=16.\n"
                   "Default kernel: B=|v-v1|/4; --strength C overrides --diameter d squared/4.\n"
                   "Without --output, calculation JSON is written to standard output.\n"
                   "--quick is for development; it overwrites production workflow outputs.\n"
                   "See docs/REPRODUCIBILITY.md for units and production settings.\n";
            return 0;
        }
        std::string command = argv[1];
        std::map<std::string, std::string> args;
        for (int i = 2; i < argc; ++i) {
            std::string key = argv[i];
            qt::require(key.rfind("--", 0) == 0, "Expected --option");
            if (key == "--quick")
                args[key] = "true";
            else {
                qt::require(i + 1 < argc, "Missing value for " + key);
                args[key] = argv[++i];
            }
        }
        auto get = [&](std::string k, std::string fallback) {
            return args.count(k) ? args.at(k) : fallback;
        };
        const std::vector<std::string> allowed = {
            "--statistics",  "--fugacity", "--basis",    "--radial",   "--angular", "--azimuth",
            "--equilibrium", "--gamma",    "--strength", "--diameter", "--output",  "--input",
            "--root",        "--quick",    "--k"};
        for (const auto &arg : args)
            qt::require(std::find(allowed.begin(), allowed.end(), arg.first) != allowed.end(),
                        "Unknown option " + arg.first);
        qt::Options o;
        o.statistics = get("--statistics", "classical");
        o.z = std::stod(get("--fugacity", "0.1"));
        o.basis = std::stoi(get("--basis", "2"));
        o.radial = std::stoi(get("--radial", "16"));
        o.angular = std::stoi(get("--angular", "12"));
        o.azimuth = std::stoi(get("--azimuth", "16"));
        o.equilibrium = std::stoi(get("--equilibrium", "96"));
        o.gamma = std::stod(get("--gamma", "1"));
        o.diameter = std::stod(get("--diameter", "1"));
        o.strength = args.count("--strength") ? std::stod(args.at("--strength"))
                                              : o.diameter * o.diameter / 4;
        std::string root = get("--root", ".");
        qt::Json out;
        if (command == "test")
            return qt::run_tests(root);
        if (command == "figures") {
            qt::make_figures(root);
            return 0;
        }
        if (command == "audit")
            out = qt::regression_audit(root);
        else if (command == "transport")
            out = qt::solve(o);
        else if (command == "reference")
            out = qt::wu_reference(o.statistics, o.z, o.radial, o.angular);
        else if (command == "wave") {
            auto r = args.count("--input") ? qt::read_json(args.at("--input")) : qt::solve(o);
            out = qt::shear_wave(r, std::stod(get("--k", "0.2")));
        } else if (command == "reproduce")
            out = qt::reproduce(root, args.count("--quick"));
        else if (command == "hilker")
            out = qt::reproduce_hilker(root, args.count("--quick"));
        else
            throw std::invalid_argument("Unknown command: " + command);
        if (args.count("--output"))
            qt::write_json(args.at("--output"), out);
        else
            std::cout << out.dump(2) << '\n';
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
