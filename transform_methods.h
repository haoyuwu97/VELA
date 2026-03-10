#ifndef VELA_TRANSFORM_METHODS_H
#define VELA_TRANSFORM_METHODS_H

#include <string>
#include <vector>

struct TrustedWindowInfo {
    int max_index;
    int trusted_points;
    double trusted_time;
    double tail_sigma;
    TrustedWindowInfo() : max_index(0), trusted_points(0), trusted_time(0.0), tail_sigma(0.0) {}
};

TrustedWindowInfo detect_trusted_window(const std::vector<double>& times,
                                        const std::vector<double>& values,
                                        int gtcut,
                                        double tail_param);

struct TransformResult {
    std::string method;
    std::vector<double> omega;
    std::vector<double> storage;
    std::vector<double> loss;
    std::vector<double> magnitude;
    std::vector<double> storage_low;
    std::vector<double> storage_high;
    std::vector<double> loss_low;
    std::vector<double> loss_high;
    double trusted_time;
    int trusted_points;
    double tail_sigma;
    std::string note;

    TransformResult() : trusted_time(0.0), trusted_points(0), tail_sigma(0.0) {}
};


TransformResult transform_irheo(const std::vector<double>& times,
                                const std::vector<double>& values,
                                int gtcut);

TransformResult transform_schwarzl(const std::vector<double>& times,
                                   const std::vector<double>& values,
                                   int gtcut);

TransformResult transform_caft_bounds(const std::vector<double>& times,
                                      const std::vector<double>& values,
                                      double dt,
                                      int raw_time_count,
                                      int gtcut,
                                      double tail_param);

#endif
