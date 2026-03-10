#include "transform_methods.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace {

struct XYData {
    std::vector<double> t;
    std::vector<double> g;
};

XYData slice_unique_monotone(const std::vector<double>& times,
                             const std::vector<double>& values,
                             int gtcut,
                             bool keep_zero) {
    XYData out;
    int limit = static_cast<int>(times.size());
    if (gtcut > 0 && gtcut < limit) {
        limit = gtcut;
    }
    if (limit <= 0) {
        return out;
    }

    double last_t = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < limit; ++i) {
        double ti = times[i];
        double gi = values[i];
        if (!std::isfinite(ti) || !std::isfinite(gi)) {
            continue;
        }
        if ((!keep_zero && ti <= 0.0) || (keep_zero && ti < 0.0)) {
            continue;
        }
        if (!out.t.empty()) {
            if (ti <= last_t) {
                if (std::fabs(ti - last_t) < 1e-14 * std::max(1.0, std::fabs(last_t))) {
                    out.g.back() = gi;
                    continue;
                }
                continue;
            }
        }
        out.t.push_back(ti);
        out.g.push_back(gi);
        last_t = ti;
    }
    return out;
}

std::vector<double> build_logspace(double start, double end, int n) {
    std::vector<double> x;
    if (!(start > 0.0) || !(end > 0.0) || n <= 0) {
        return x;
    }
    if (n == 1) {
        x.push_back(start);
        return x;
    }
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    x.resize(n);
    for (int i = 0; i < n; ++i) {
        double s = static_cast<double>(i) / static_cast<double>(n - 1);
        x[i] = std::pow(10.0, log_start + (log_end - log_start) * s);
    }
    return x;
}

std::vector<double> build_default_omega_grid(const std::vector<double>& times) {
    std::vector<double> omega;
    if (times.size() < 2) {
        return omega;
    }
    for (std::size_t i = times.size() - 1; i >= 1; --i) {
        if (times[i] > 0.0) {
            omega.push_back(1.0 / times[i]);
        }
        if (i == 1) {
            break;
        }
    }
    return omega;
}

void piecewise_linear_transform(const std::vector<double>& times,
                                const std::vector<double>& values,
                                const std::vector<double>& omega,
                                std::vector<double>& storage,
                                std::vector<double>& loss) {
    storage.assign(omega.size(), 0.0);
    loss.assign(omega.size(), 0.0);
    if (times.size() < 2 || values.size() != times.size()) {
        return;
    }
    for (std::size_t iw = 0; iw < omega.size(); ++iw) {
        double w = omega[iw];
        if (!(w > 0.0)) {
            continue;
        }
        double w2 = w * w;
        double gp = 0.0;
        double gpp = 0.0;
        for (std::size_t j = 1; j < times.size(); ++j) {
            double t0 = times[j - 1];
            double t1 = times[j];
            double g0 = values[j - 1];
            double g1 = values[j];
            double dtk = t1 - t0;
            if (!(dtk > 0.0)) {
                continue;
            }
            double cost0 = std::cos(w * t0);
            double sint0 = std::sin(w * t0);
            double coswdtk = std::cos(w * dtk);
            double sinwdtk = std::sin(w * dtk);
            double X1 = g1 * (((coswdtk - 1.0) / (w2 * dtk)) + sinwdtk / w)
                      - g0 * ((coswdtk - 1.0) / (w2 * dtk));
            double X2 = g1 * (coswdtk / w - sinwdtk / (w2 * dtk))
                      - g0 * (1.0 / w - sinwdtk / (w2 * dtk));
            gp += sint0 * X1 - cost0 * X2;
            gpp += sint0 * X2 + cost0 * X1;
        }
        storage[iw] = gp * w;
        loss[iw] = gpp * w;
    }
}

void fill_magnitude_and_bounds(TransformResult& result) {
    result.magnitude.resize(result.omega.size(), 0.0);
    if (result.storage_low.empty()) {
        result.storage_low = result.storage;
    }
    if (result.storage_high.empty()) {
        result.storage_high = result.storage;
    }
    if (result.loss_low.empty()) {
        result.loss_low = result.loss;
    }
    if (result.loss_high.empty()) {
        result.loss_high = result.loss;
    }
    for (std::size_t i = 0; i < result.omega.size(); ++i) {
        result.magnitude[i] = std::sqrt(result.storage[i] * result.storage[i] + result.loss[i] * result.loss[i]);
        if (result.storage_low[i] > result.storage_high[i]) {
            std::swap(result.storage_low[i], result.storage_high[i]);
        }
        if (result.loss_low[i] > result.loss_high[i]) {
            std::swap(result.loss_low[i], result.loss_high[i]);
        }
    }
}

std::vector<double> interpolate_linear(const std::vector<double>& x_in,
                                       const std::vector<double>& y_in,
                                       const std::vector<double>& x_out) {
    std::vector<double> y_out(x_out.size(), 0.0);
    if (x_in.empty()) {
        return y_out;
    }
    if (x_in.size() == 1) {
        std::fill(y_out.begin(), y_out.end(), y_in[0]);
        return y_out;
    }
    std::size_t j = 0;
    for (std::size_t i = 0; i < x_out.size(); ++i) {
        double x = x_out[i];
        while (j + 1 < x_in.size() && x_in[j + 1] < x) {
            ++j;
        }
        if (x <= x_in.front()) {
            double dx = x_in[1] - x_in[0];
            double slope = (y_in[1] - y_in[0]) / dx;
            y_out[i] = y_in[0] + slope * (x - x_in[0]);
        } else if (x >= x_in.back()) {
            std::size_t n = x_in.size();
            double dx = x_in[n - 1] - x_in[n - 2];
            double slope = (y_in[n - 1] - y_in[n - 2]) / dx;
            y_out[i] = y_in[n - 1] + slope * (x - x_in[n - 1]);
        } else {
            double dx = x_in[j + 1] - x_in[j];
            double alpha = (x - x_in[j]) / dx;
            y_out[i] = (1.0 - alpha) * y_in[j] + alpha * y_in[j + 1];
        }
    }
    return y_out;
}

std::vector<double> polynomial_weights_at_zero(const std::vector<double>& x) {
    std::vector<double> weights(x.size(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        double wi = 1.0;
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (j == i) {
                continue;
            }
            wi *= (0.0 - x[j]) / (x[i] - x[j]);
        }
        weights[i] = wi;
    }
    return weights;
}

double extrapolate_g0(const std::vector<double>& t, const std::vector<double>& g) {
    if (t.empty()) {
        return 0.0;
    }
    if (std::fabs(t.front()) < 1e-14) {
        return g.front();
    }
    std::size_t order = std::min<std::size_t>(4, t.size());
    std::vector<double> tx(order), gy(order);
    for (std::size_t i = 0; i < order; ++i) {
        tx[i] = t[i];
        gy[i] = g[i];
    }
    std::vector<double> w = polynomial_weights_at_zero(tx);
    double g0 = 0.0;
    for (std::size_t i = 0; i < order; ++i) {
        g0 += w[i] * gy[i];
    }
    return g0;
}

void build_range_with_zero(const XYData& data,
                           std::vector<double>& t,
                           std::vector<double>& g) {
    t = data.t;
    g = data.g;
    if (t.empty()) {
        return;
    }
    if (t.front() > 0.0) {
        double g0 = extrapolate_g0(t, g);
        t.insert(t.begin(), 0.0);
        g.insert(g.begin(), g0);
    }
}

std::vector<double> matvec_exp(const std::vector<double>& time_shift,
                               const std::vector<double>& tau,
                               const std::vector<double>& amp) {
    std::vector<double> y(time_shift.size(), 0.0);
    for (std::size_t i = 0; i < time_shift.size(); ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k < tau.size(); ++k) {
            s += amp[k] * std::exp(-time_shift[i] / tau[k]);
        }
        y[i] = s;
    }
    return y;
}

std::vector<double> fit_nonnegative_tail(const std::vector<double>& time_shift,
                                         const std::vector<double>& y_target,
                                         const std::vector<double>& tau,
                                         double lambda) {
    std::size_t m = time_shift.size();
    std::size_t n = tau.size();
    std::vector<double> amp(n, 0.0);
    if (m == 0 || n == 0) {
        return amp;
    }

    std::vector< std::vector<double> > A(m, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t k = 0; k < n; ++k) {
            A[i][k] = std::exp(-time_shift[i] / tau[k]);
        }
    }

    std::vector<double> v(n, 1.0 / std::sqrt(static_cast<double>(n)));
    for (int it = 0; it < 20; ++it) {
        std::vector<double> Av(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t k = 0; k < n; ++k) {
                Av[i] += A[i][k] * v[k];
            }
        }
        std::vector<double> AtAv(n, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t i = 0; i < m; ++i) {
                AtAv[k] += A[i][k] * Av[i];
            }
            AtAv[k] += lambda * v[k];
        }
        double norm = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            norm += AtAv[k] * AtAv[k];
        }
        norm = std::sqrt(std::max(norm, 1e-30));
        for (std::size_t k = 0; k < n; ++k) {
            v[k] = AtAv[k] / norm;
        }
    }
    std::vector<double> Av(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t k = 0; k < n; ++k) {
            Av[i] += A[i][k] * v[k];
        }
    }
    std::vector<double> AtAv(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t i = 0; i < m; ++i) {
            AtAv[k] += A[i][k] * Av[i];
        }
        AtAv[k] += lambda * v[k];
    }
    double lipschitz = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        lipschitz += v[k] * AtAv[k];
    }
    lipschitz = std::max(lipschitz, 1e-8);
    double alpha = 0.9 / lipschitz;

    for (int iter = 0; iter < 5000; ++iter) {
        std::vector<double> model(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t k = 0; k < n; ++k) {
                model[i] += A[i][k] * amp[k];
            }
        }
        std::vector<double> grad(n, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t i = 0; i < m; ++i) {
                grad[k] += A[i][k] * (model[i] - y_target[i]);
            }
            grad[k] += lambda * amp[k];
        }
        double max_delta = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            double old = amp[k];
            amp[k] -= alpha * grad[k];
            if (amp[k] < 0.0) {
                amp[k] = 0.0;
            }
            max_delta = std::max(max_delta, std::fabs(amp[k] - old));
        }
        if (max_delta < 1e-12) {
            break;
        }
    }
    return amp;
}

TransformResult package_piecewise_result(const std::string& name,
                                         const std::vector<double>& times,
                                         const std::vector<double>& values,
                                         const std::vector<double>& omega,
                                         const std::string& note) {
    TransformResult out;
    out.method = name;
    out.omega = omega;
    piecewise_linear_transform(times, values, omega, out.storage, out.loss);
    out.note = note;
    fill_magnitude_and_bounds(out);
    return out;
}

} // namespace

TrustedWindowInfo detect_trusted_window(const std::vector<double>& times,
                                        const std::vector<double>& values,
                                        int gtcut,
                                        double tail_param) {
    TrustedWindowInfo info;
    int limit = static_cast<int>(values.size());
    if (gtcut > 0 && gtcut < limit) {
        limit = gtcut;
    }
    if (limit <= 0) {
        return info;
    }
    info.max_index = limit;
    if (limit <= 2) {
        info.trusted_points = limit;
        info.trusted_time = times[limit - 1];
        return info;
    }
    int window = std::max(5, limit / 10);
    int tail_start = std::max(0, limit - window);
    double mean = 0.0;
    for (int i = tail_start; i < limit; ++i) {
        mean += values[i];
    }
    mean /= static_cast<double>(limit - tail_start);
    double var = 0.0;
    for (int i = tail_start; i < limit; ++i) {
        double diff = values[i] - mean;
        var += diff * diff;
    }
    var /= static_cast<double>(std::max(1, limit - tail_start - 1));
    info.tail_sigma = std::sqrt(std::max(var, 0.0));
    double k = (tail_param > 0.0 && std::isfinite(tail_param)) ? tail_param : 2.0;
    int consecutive = std::max(3, window / 4);
    int run = 0;
    int cutoff = limit;
    for (int i = limit - 1; i >= 0; --i) {
        if (std::fabs(values[i]) <= k * info.tail_sigma) {
            ++run;
            if (run >= consecutive) {
                cutoff = i;
                break;
            }
        } else {
            run = 0;
        }
    }
    if (cutoff < limit) {
        limit = std::max(cutoff, 2);
    }
    info.max_index = limit;
    info.trusted_points = limit;
    info.trusted_time = times[limit - 1];
    return info;
}

TransformResult transform_irheo(const std::vector<double>& times,
                                const std::vector<double>& values,
                                int gtcut) {
    TransformResult out;
    out.method = "irheo";

    XYData data = slice_unique_monotone(times, values, gtcut, true);
    if (data.t.size() < 2) {
        out.note = "Too few unique data points for i-RheoFT";
        return out;
    }

    double g0 = extrapolate_g0(data.t, data.g);
    std::size_t ind1 = 0;
    while (ind1 < data.t.size() && !(data.t[ind1] > 0.0)) {
        ++ind1;
    }
    if (ind1 >= data.t.size()) {
        out.note = "No positive times available for i-RheoFT";
        return out;
    }

    double t1 = data.t[ind1];
    double g1 = data.g[ind1];
    double tinf = data.t.back();
    if (!(t1 > 0.0) || !(tinf > 0.0)) {
        out.note = "Invalid positive time range for i-RheoFT";
        return out;
    }

    int nfreq = static_cast<int>(data.t.size());
    out.omega = build_logspace(1.0 / tinf, 1.0 / t1, nfreq);
    out.storage.assign(out.omega.size(), 0.0);
    out.loss.assign(out.omega.size(), 0.0);

    std::vector<double> coeff;
    if (data.t.size() >= ind1 + 2) {
        coeff.resize(data.t.size() - ind1 - 1, 0.0);
        for (std::size_t i = ind1; i + 1 < data.t.size(); ++i) {
            coeff[i - ind1] = (data.g[i + 1] - data.g[i]) / (data.t[i + 1] - data.t[i]);
        }
    }

    for (std::size_t iw = 0; iw < out.omega.size(); ++iw) {
        double w = out.omega[iw];
        if (!(w > 0.0)) {
            continue;
        }
        double gp = g0 + std::sin(w * t1) * (g1 - g0) / (w * t1);
        double gpp = -(1.0 - std::cos(w * t1)) * (g1 - g0) / (w * t1);
        for (std::size_t i = ind1; i + 1 < data.t.size(); ++i) {
            double slope = coeff[i - ind1];
            gp += slope * (-std::sin(w * data.t[i]) + std::sin(w * data.t[i + 1])) / w;
            gpp += -slope * (std::cos(w * data.t[i]) - std::cos(w * data.t[i + 1])) / w;
        }
        out.storage[iw] = gp;
        out.loss[iw] = gpp;
    }

    out.trusted_points = static_cast<int>(data.t.size());
    out.trusted_time = tinf;
    out.note = "i-RheoFT baseline ported from the open RepTate G(t) application (non-oversampled viewiRheo logic).";
    fill_magnitude_and_bounds(out);
    return out;
}

TransformResult transform_schwarzl(const std::vector<double>& times,
                                   const std::vector<double>& values,
                                   int gtcut) {
    TransformResult out;
    out.method = "schwarzl";

    XYData data0 = slice_unique_monotone(times, values, gtcut, true);
    std::vector<double> t;
    std::vector<double> g;
    build_range_with_zero(data0, t, g);
    if (t.size() < 2) {
        out.note = "Too few data points for Schwarzl baseline";
        return out;
    }

    double t1 = 0.0;
    std::size_t first_pos = 0;
    while (first_pos < t.size() && !(t[first_pos] > 0.0)) {
        ++first_pos;
    }
    if (first_pos >= t.size()) {
        out.note = "No positive times available for Schwarzl baseline";
        return out;
    }
    t1 = t[first_pos];
    double tinf = t.back();
    int nfreq = std::max(16, static_cast<int>(data0.t.size()));
    out.omega = build_logspace(1.0 / tinf, 1.0 / t1, nfreq);

    int nresamp = std::max(64, 4 * static_cast<int>(data0.t.size()));
    std::vector<double> xgrid = build_logspace(t1, tinf, nresamp);
    std::vector<double> ygrid = interpolate_linear(data0.t, data0.g, xgrid);
    xgrid.insert(xgrid.begin(), 0.0);
    ygrid.insert(ygrid.begin(), extrapolate_g0(t, g));

    piecewise_linear_transform(xgrid, ygrid, out.omega, out.storage, out.loss);
    out.trusted_points = static_cast<int>(data0.t.size());
    out.trusted_time = tinf;
    out.note = "Schwarzl-style numerical baseline on a log-resampled time grid. This is a transparent numerical baseline, not a byte-identical reproduction of RepTate's helper library.";
    fill_magnitude_and_bounds(out);
    return out;
}

TransformResult transform_caft_bounds(const std::vector<double>& times,
                                      const std::vector<double>& values,
                                      double dt,
                                      int /*raw_time_count*/,
                                      int gtcut,
                                      double tail_param) {
    TransformResult out;
    out.method = "caft";

    TrustedWindowInfo info = detect_trusted_window(times, values, gtcut, tail_param);
    if (info.trusted_points < 2) {
        out.note = "Too few trusted points for CAFT-Bounds";
        return out;
    }

    int trusted = info.trusted_points;
    out.trusted_points = trusted;
    out.trusted_time = info.trusted_time;
    out.tail_sigma = info.tail_sigma;

    int fit_count = std::max(6, trusted / 3);
    fit_count = std::min(fit_count, trusted);
    int fit_start = trusted - fit_count;
    if (fit_start < 0) {
        fit_start = 0;
    }

    std::vector<double> fit_t;
    std::vector<double> fit_y_center;
    std::vector<double> fit_y_low;
    std::vector<double> fit_y_high;
    fit_t.reserve(fit_count);
    fit_y_center.reserve(fit_count);
    fit_y_low.reserve(fit_count);
    fit_y_high.reserve(fit_count);

    for (int i = fit_start; i < trusted; ++i) {
        fit_t.push_back(times[i] - times[fit_start]);
        double y = values[i];
        fit_y_center.push_back(std::max(y, 0.0));
        fit_y_low.push_back(std::max(y - info.tail_sigma, 0.0));
        fit_y_high.push_back(std::max(y + info.tail_sigma, 0.0));
    }

    double t_end = times[trusted - 1];
    double t_last = times[std::min<int>(static_cast<int>(times.size()), gtcut > 0 ? gtcut : static_cast<int>(times.size())) - 1];
    double tau_min = dt;
    if (fit_t.size() >= 2) {
        tau_min = std::max(dt, 0.5 * (fit_t[1] - fit_t[0]));
    }
    double tau_max = std::max(t_last, t_end) * 20.0;
    if (!(tau_max > tau_min)) {
        tau_max = tau_min * 10.0;
    }
    const int K = 16;
    std::vector<double> tau = build_logspace(tau_min, tau_max, K);
    double y_scale = 0.0;
    for (std::size_t i = 0; i < fit_y_center.size(); ++i) {
        y_scale = std::max(y_scale, fit_y_center[i]);
    }
    double lambda = std::max(1e-12, 1e-6 * std::max(1.0, y_scale * y_scale));

    std::vector<double> amp_center = fit_nonnegative_tail(fit_t, fit_y_center, tau, lambda);
    std::vector<double> amp_low = fit_nonnegative_tail(fit_t, fit_y_low, tau, lambda);
    std::vector<double> amp_high = fit_nonnegative_tail(fit_t, fit_y_high, tau, lambda);

    double target_tc = std::max(values[trusted - 1], 0.0);
    double tshift_tc = times[trusted - 1] - times[fit_start];
    double model_tc_center = 0.0;
    double model_tc_low = 0.0;
    double model_tc_high = 0.0;
    for (int k = 0; k < K; ++k) {
        double fac = std::exp(-tshift_tc / tau[k]);
        model_tc_center += amp_center[k] * fac;
        model_tc_low += amp_low[k] * fac;
        model_tc_high += amp_high[k] * fac;
    }
    if (model_tc_center > 0.0 && target_tc > 0.0) {
        double s = target_tc / model_tc_center;
        for (int k = 0; k < K; ++k) amp_center[k] *= s;
    }
    if (model_tc_low > 0.0 && target_tc > 0.0) {
        double s = std::max(target_tc - info.tail_sigma, 0.0) / model_tc_low;
        for (int k = 0; k < K; ++k) amp_low[k] *= s;
    }
    if (model_tc_high > 0.0 && target_tc > 0.0) {
        double s = (target_tc + info.tail_sigma) / model_tc_high;
        for (int k = 0; k < K; ++k) amp_high[k] *= s;
    }

    double tail_end = std::max(t_last, tau_max * 6.0);
    int nsynth = 96;
    double start_tail = std::max(t_end * (1.0 + 1e-6), t_end + dt);
    std::vector<double> tail_times = build_logspace(start_tail, tail_end, nsynth);
    std::vector<double> tail_shift(tail_times.size(), 0.0);
    for (std::size_t i = 0; i < tail_times.size(); ++i) {
        tail_shift[i] = tail_times[i] - times[fit_start];
    }
    std::vector<double> tail_center = matvec_exp(tail_shift, tau, amp_center);
    std::vector<double> tail_low = matvec_exp(tail_shift, tau, amp_low);
    std::vector<double> tail_high = matvec_exp(tail_shift, tau, amp_high);

    std::vector<double> t_comb_center(times.begin(), times.begin() + trusted);
    std::vector<double> g_comb_center(values.begin(), values.begin() + trusted);
    std::vector<double> t_comb_low = t_comb_center;
    std::vector<double> g_comb_low = g_comb_center;
    std::vector<double> t_comb_high = t_comb_center;
    std::vector<double> g_comb_high = g_comb_center;

    for (std::size_t i = 0; i < tail_times.size(); ++i) {
        t_comb_center.push_back(tail_times[i]);
        g_comb_center.push_back(tail_center[i]);
        t_comb_low.push_back(tail_times[i]);
        g_comb_low.push_back(tail_low[i]);
        t_comb_high.push_back(tail_times[i]);
        g_comb_high.push_back(tail_high[i]);
    }

    out.omega = build_logspace(1.0 / tail_end, 1.0 / std::max(dt, times[1]), std::max(32, trusted));
    piecewise_linear_transform(t_comb_center, g_comb_center, out.omega, out.storage, out.loss);

    std::vector<double> gp_low_a, gpp_low_a, gp_high_a, gpp_high_a;
    piecewise_linear_transform(t_comb_low, g_comb_low, out.omega, gp_low_a, gpp_low_a);
    piecewise_linear_transform(t_comb_high, g_comb_high, out.omega, gp_high_a, gpp_high_a);

    out.storage_low.resize(out.omega.size(), 0.0);
    out.storage_high.resize(out.omega.size(), 0.0);
    out.loss_low.resize(out.omega.size(), 0.0);
    out.loss_high.resize(out.omega.size(), 0.0);
    for (std::size_t i = 0; i < out.omega.size(); ++i) {
        out.storage_low[i] = std::min(gp_low_a[i], gp_high_a[i]);
        out.storage_high[i] = std::max(gp_low_a[i], gp_high_a[i]);
        out.loss_low[i] = std::min(gpp_low_a[i], gpp_high_a[i]);
        out.loss_high[i] = std::max(gpp_low_a[i], gpp_high_a[i]);
        out.storage_low[i] = std::min(out.storage_low[i], out.storage[i]);
        out.storage_high[i] = std::max(out.storage_high[i], out.storage[i]);
        out.loss_low[i] = std::min(out.loss_low[i], out.loss[i]);
        out.loss_high[i] = std::max(out.loss_high[i], out.loss[i]);
    }

    std::ostringstream oss;
    oss << "CAFT-Bounds v1: trusted-prefix detection plus nonnegative multi-exponential tail continuation; bounds are generated from ±sigma tail envelopes on the trusted window.";
    out.note = oss.str();
    fill_magnitude_and_bounds(out);
    return out;
}
