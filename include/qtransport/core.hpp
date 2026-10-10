#pragma once
#include <Eigen/Dense>
#include <array>
#include <complex>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace qt {
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using CMat = Eigen::MatrixXcd;
using CVec = Eigen::VectorXcd;
using Json = nlohmann::json;
constexpr double pi = 3.1415926535897932384626433832795;
constexpr double kb = 1.380649e-23, hbar = 1.0545718176461565e-34;
constexpr double a0 = 5.29177210903e-11, mass_k39 = 38.9637064864 * 1.66053906660e-27;
void require(bool condition, const std::string &message);
Json array(const Vec &v);
Json array(const Mat &m);
Vec vector(const Json &j);
Mat matrix(const Json &j);
Vec linspace(double a, double b, int n);
Json read_json(const std::string &path);
void write_json(const std::string &path, const Json &value);
struct Rule {
    Vec x, w;
};
Rule laguerre_rule(int n, double alpha);
Rule legendre_rule(int n);
Vec laguerre_values(double x, int count, double alpha);
struct Moments {
    Vec shifts, dh, ds, matching;
    Mat hh, hs;
};
struct Equilibrium {
    std::string statistics;
    double z, g3_over_z = 0, normalization = 0, occupation_conversion = 0, pressure = 0;
    int eps, order;
    Vec x, f_measure, w_measure;
    Equilibrium(std::string statistics, double z, int order = 96);
    double denominator(double energy) const;
    double density(const Eigen::Vector3d &v) const;
    Vec shifts(int count) const;
    std::pair<Vec, Vec> profiles(double energy, int count, const Vec &shift) const;
    Moments moments(int count) const;
};
struct Options {
    std::string statistics = "classical";
    double z = .1, diameter = 1, gamma = 1, strength = .25;
    int basis = 2, radial = 16, angular = 12, azimuth = 16, equilibrium = 96;
};
struct Brackets {
    Mat heat, shear;
    Json diagnostics;
};
Brackets brackets(const Equilibrium &eq, const Options &o);
Json relaxation(const Mat &h, const Mat &m, const Vec &d, const Vec &times);
Json solve(const Options &o);
Json wu_reference(const std::string &statistics, double z, int radial = 32, int angular = 24);
Json shear_wave(const Json &result, double wave_number = .2, const Vec &times = Vec());
std::array<double, 3> bose_functions(double z);
Json normal_state(double n, double temperature, double a, double mass = mass_k39);
std::pair<double, double> dimensional_transport(const Json &r, double temperature, double a,
                                                double mass = mass_k39);
CVec sound_modes(double qe, double qk, double ratio);
Json sound_prediction(const Json &state, double eta, double kappa, double length);
CMat longitudinal_matrix(const Json &state, double eta, double kappa, double k);
Json observations(const std::string &csv_path);
Json read_table(const std::string &csv_path);
Json state_from_observation(const Json &row, bool mean_field = true);
Json representation(const Json &result, double half_width, int points);
Json reproduce(const std::string &root, bool quick);
Json reproduce_hilker(const std::string &root, bool quick);
void make_figures(const std::string &root);
Json regression_audit(const std::string &root);
int run_tests(const std::string &root);
} // namespace qt
