#include "qtransport/core.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>

namespace qt {
void require(bool ok, const std::string &message) {
    if (!ok)
        throw std::invalid_argument(message);
}
Json array(const Vec &v) {
    Json a = Json::array();
    for (Eigen::Index i = 0; i < v.size(); ++i)
        a.push_back(v[i]);
    return a;
}
Json array(const Mat &m) {
    Json a = Json::array();
    for (Eigen::Index i = 0; i < m.rows(); ++i) {
        Json row = Json::array();
        for (Eigen::Index j = 0; j < m.cols(); ++j)
            row.push_back(m(i, j));
        a.push_back(row);
    }
    return a;
}
Vec vector(const Json &j) {
    Vec v(j.size());
    for (size_t i = 0; i < j.size(); ++i)
        v[i] = j[i].get<double>();
    return v;
}
Mat matrix(const Json &j) {
    require(j.is_array() && !j.empty(), "Empty matrix");
    Mat m(j.size(), j[0].size());
    for (size_t i = 0; i < j.size(); ++i) {
        require(j[i].size() == size_t(m.cols()), "Ragged matrix");
        for (size_t k = 0; k < j[i].size(); ++k)
            m(i, k) = j[i][k].get<double>();
    }
    return m;
}
Vec linspace(double a, double b, int n) {
    require(n >= 2, "At least two time points required");
    return Vec::LinSpaced(n, a, b);
}
Json read_json(const std::string &path) {
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Cannot read " + path);
    Json j;
    f >> j;
    return j;
}
void write_json(const std::string &path, const Json &j) {
    auto p = std::filesystem::path(path);
    if (p.has_parent_path())
        std::filesystem::create_directories(p.parent_path());
    std::ofstream f(path);
    if (!f)
        throw std::runtime_error("Cannot write " + path);
    f << j.dump(2) << '\n';
}
Rule laguerre_rule(int n, double alpha) {
    require(n >= 1 && std::isfinite(alpha) && alpha > -1, "Invalid Laguerre rule");
    static std::map<std::pair<int, double>, Rule> cache;
    auto key = std::make_pair(n, alpha);
    if (cache.count(key))
        return cache.at(key);
    Vec diag(n), sub(n - 1);
    for (int k = 0; k < n; ++k)
        diag[k] = 2 * k + 1 + alpha;
    for (int k = 1; k < n; ++k)
        sub[k - 1] = std::sqrt(k * (k + alpha));
    Eigen::SelfAdjointEigenSolver<Mat> eig;
    eig.computeFromTridiagonal(diag, sub);
    if (eig.info() != Eigen::Success)
        throw std::runtime_error("Laguerre eigensolver failed");
    Rule r{eig.eigenvalues(), std::tgamma(alpha + 1) *
                                  eig.eigenvectors().row(0).array().square().matrix().transpose()};
    cache.emplace(key, r);
    return r;
}
Rule legendre_rule(int n) {
    require(n >= 1, "Invalid Legendre rule");
    Vec diag = Vec::Zero(n), sub(n - 1);
    for (int k = 1; k < n; ++k)
        sub[k - 1] = k / std::sqrt(4.0 * k * k - 1);
    Eigen::SelfAdjointEigenSolver<Mat> eig;
    eig.computeFromTridiagonal(diag, sub);
    if (eig.info() != Eigen::Success)
        throw std::runtime_error("Legendre eigensolver failed");
    return {eig.eigenvalues(), 2 * eig.eigenvectors().row(0).array().square().matrix().transpose()};
}
Vec laguerre_values(double x, int count, double alpha) {
    require(count > 0, "Empty polynomial basis");
    Vec p(count);
    p[0] = 1;
    if (count > 1)
        p[1] = 1 + alpha - x;
    for (int n = 1; n < count - 1; ++n)
        p[n + 1] = ((2 * n + 1 + alpha - x) * p[n] - (n + alpha) * p[n - 1]) / (n + 1);
    return p;
}
Equilibrium::Equilibrium(std::string stat, double fugacity, int radial)
    : statistics(std::move(stat)), z(fugacity), order(radial) {
    require(statistics == "classical" || statistics == "bose" || statistics == "fermi",
            "Unknown statistics");
    eps = statistics == "bose" ? 1 : (statistics == "fermi" ? -1 : 0);
    require(std::isfinite(z) && z > 0 && (eps != 1 || z < 1), "Require z>0 and normal Bose z<1");
    auto rule = laguerre_rule(order, .5);
    x = rule.x;
    Vec den(order);
    for (int i = 0; i < order; ++i)
        den[i] = denominator(x[i]);
    Vec maxwell = rule.w / std::tgamma(1.5);
    g3_over_z = (maxwell.array() / den.array()).sum();
    normalization = std::pow(2 * pi, 1.5) * g3_over_z;
    occupation_conversion = z * normalization;
    f_measure = (maxwell.array() / (g3_over_z * den.array())).matrix();
    w_measure = (f_measure.array() / den.array()).matrix();
    pressure = f_measure.dot(2 * x) / 3;
}
double Equilibrium::denominator(double e) const {
    return 1 - eps * z * std::exp(-e);
}
double Equilibrium::density(const Eigen::Vector3d &v) const {
    double e = v.squaredNorm() / 2;
    return std::exp(-e) / (normalization * denominator(e));
}
Vec Equilibrium::shifts(int count) const {
    Vec out = Vec::Zero(count);
    double total = 0;
    for (int i = 0; i < order; ++i) {
        double w = w_measure[i] * 2 * x[i];
        Vec p = laguerre_values(x[i], count + 1, 1.5);
        out += w * p.tail(count);
        total += w;
    }
    return out / total;
}
std::pair<Vec, Vec> Equilibrium::profiles(double e, int count, const Vec &shift) const {
    return {laguerre_values(e, count + 1, 1.5).tail(count) - shift, laguerre_values(e, count, 2.5)};
}
Moments Equilibrium::moments(int count) const {
    Moments m{shifts(count),    Vec::Zero(count),        Vec::Zero(count),
              Vec::Zero(count), Mat::Zero(count, count), Mat::Zero(count, count)};
    for (int i = 0; i < order; ++i) {
        auto [p, q] = profiles(x[i], count, m.shifts);
        double v2 = 2 * x[i], wh = w_measure[i] * v2 / 3, ws = w_measure[i] * v2 * v2 / 15;
        m.dh += wh * (x[i] - 2.5 * pressure) * p;
        m.ds += ws * q;
        m.hh.noalias() += wh * (p * p.transpose());
        m.hs.noalias() += ws * (q * q.transpose());
        m.matching += 3 * wh * p;
    }
    return m;
}
} // namespace qt
