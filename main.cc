#include "correlator.h"
#include "transform_methods.h"
#include "lammps_mua_backend.h"

#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct MuAStats {
    bool provided;
    bool from_series;
    double mean_raw;
    double stderr_raw;
    int count;
    std::string source;
    MuAStats() : provided(false), from_series(false), mean_raw(0.0), stderr_raw(0.0), count(0), source("none") {}
};

void print_usage() {
    std::cout << "Usage:" << std::endl
              << "  VELA INPUT_FILE OUTPUT_FILE_Gt OUTPUT_FILE_Gw DT MODE [GTCUT] [TAIL_PARAM] [--muA VALUE] [--muA-series FILE] [--muA-lammps-input FILE] [--muA-lammps-lib PATH] [--muA-lammps-numdiff DELTA] [--muA-lammps-channel CHANNEL] [--transform METHOD]" << std::endl
              << "Notes:" << std::endl
              << "  VELA outputs raw reduced/absolute moduli. For LJ/reduced units use V/kBT; for LAMMPS pressure-unit input use V/(kB*T*nktv2p)." << std::endl
              << "  MODE=0 uses input columns xy,xz,yz." << std::endl
              << "  MODE=1 uses input columns xx,yy,zz,xy,xz,yz and internally forms normal-stress differences plus shear stresses." << std::endl
              << "  GTCUT sets the maximum time index used for transform-range selection or trusted-window detection; the default is the full time range." << std::endl
              << "  TAIL_PARAM is the k-sigma threshold used by the trusted-window detector." << std::endl
              << "  --muA VALUE supplies a raw affine modulus and switches to the absolute branch." << std::endl
              << "  --muA-series FILE supplies a raw muA time/block series and switches to the absolute branch." << std::endl
              << "  --muA-lammps-input FILE asks VELA to call liblammps and compute muA from a complete LAMMPS setup script." << std::endl
              << "  --muA-lammps-lib PATH optionally specifies the shared liblammps path." << std::endl
              << "  --muA-lammps-numdiff DELTA switches the born/matrix backend to numdiff with the supplied strain delta." << std::endl
              << "  --muA-lammps-channel CHANNEL chooses auto|xy|xz|yz|avg3|iso6 (default: auto)." << std::endl
              << "  --transform METHOD chooses direct|irheo|schwarzl|caft (default: direct)." << std::endl;
}

bool parse_double(const std::string &s, double &value) {
    char *end = NULL;
    value = std::strtod(s.c_str(), &end);
    return end != s.c_str() && end != NULL && *end == '\0' && std::isfinite(value);
}

bool parse_int(const std::string &s, int &value) {
    char *end = NULL;
    long tmp = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || end == NULL || *end != '\0') {
        return false;
    }
    value = static_cast<int>(tmp);
    return true;
}

std::string normalize_transform(std::string method) {
    for (std::size_t i = 0; i < method.size(); ++i) {
        method[i] = static_cast<char>(std::tolower(method[i]));
    }
    if (method == "i-rheo" || method == "irheoft" || method == "i-rheoft") {
        return "irheo";
    }
    if (method == "schwarzl-baseline" || method == "schwarzl_baseline") {
        return "schwarzl";
    }
    if (method == "caft-bounds" || method == "caft_bounds") {
        return "caft";
    }
    return method;
}

bool load_muA_series(const std::string& path, MuAStats& stats) {
    std::ifstream fin(path.c_str());
    if (!fin.is_open() || fin.bad()) {
        return false;
    }
    std::string line;
    std::vector<double> vals;
    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }
        std::size_t hash_pos = line.find('#');
        if (hash_pos != std::string::npos) {
            line = line.substr(0, hash_pos);
        }
        std::stringstream ss(line);
        std::vector<double> nums;
        std::string tok;
        while (ss >> tok) {
            double x = 0.0;
            if (parse_double(tok, x)) {
                nums.push_back(x);
            }
        }
        if (!nums.empty()) {
            vals.push_back(nums.back());
        }
    }
    if (vals.empty()) {
        return false;
    }
    double mean = 0.0;
    for (std::size_t i = 0; i < vals.size(); ++i) {
        mean += vals[i];
    }
    mean /= static_cast<double>(vals.size());
    double var = 0.0;
    for (std::size_t i = 0; i < vals.size(); ++i) {
        double diff = vals[i] - mean;
        var += diff * diff;
    }
    if (vals.size() > 1) {
        var /= static_cast<double>(vals.size() - 1);
    } else {
        var = 0.0;
    }
    stats.provided = true;
    stats.from_series = true;
    stats.mean_raw = mean;
    stats.stderr_raw = (vals.size() > 1) ? std::sqrt(var / static_cast<double>(vals.size())) : 0.0;
    stats.count = static_cast<int>(vals.size());
    stats.source = path;
    return true;
}

void apply_absolute_shift(TransformResult& tr, double Ge_raw) {
    for (std::size_t i = 0; i < tr.storage.size(); ++i) {
        tr.storage[i] += Ge_raw;
        if (!tr.storage_low.empty()) {
            tr.storage_low[i] += Ge_raw;
        }
        if (!tr.storage_high.empty()) {
            tr.storage_high[i] += Ge_raw;
        }
        tr.magnitude[i] = std::sqrt(tr.storage[i] * tr.storage[i] + tr.loss[i] * tr.loss[i]);
    }
}

TransformResult collect_direct_result(Correlator& c, double dt, const TrustedWindowInfo& info) {
    TransformResult tr;
    tr.method = "direct";
    tr.trusted_time = info.trusted_time;
    tr.trusted_points = info.trusted_points;
    tr.tail_sigma = info.tail_sigma;
    tr.note = "Current VELA direct transform on the multi-tau grid with the legacy k-sigma confidence cutoff reinterpreted as a trusted-window detector.";
    if (c.npcorr <= 1) {
        return tr;
    }
    for (int i = static_cast<int>(c.npcorr) - 1; i > 0; --i) {
        double omega = 1.0 / (c.wt[i] * dt);
        tr.omega.push_back(omega);
        tr.storage.push_back(c.Gw_storage[i]);
        tr.loss.push_back(c.Gw_loss[i]);
        tr.magnitude.push_back(c.Gw[i]);
        tr.storage_low.push_back(c.Gw_storage[i]);
        tr.storage_high.push_back(c.Gw_storage[i]);
        tr.loss_low.push_back(c.Gw_loss[i]);
        tr.loss_high.push_back(c.Gw_loss[i]);
    }
    return tr;
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc < 6) {
        print_usage();
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_gt_file = argv[2];
    const char *output_gw_file = argv[3];

    double dt = std::atof(argv[4]);
    int mode = std::atoi(argv[5]);
    int Gtcut = -1;
    double tail_param = 2.0;
    std::string transform_method = "direct";
    MuAStats muA_stats;

    bool have_gtcut = false;
    bool have_tail = false;
    bool muA_scalar_seen = false;
    bool muA_series_seen = false;
    bool muA_lammps_seen = false;
    LammpsMuAOptions muA_lammps_options;
    LammpsMuAResult muA_lammps_result;

    for (int i = 6; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--muA") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --muA REQUIRES A VALUE!" << std::endl;
                return 1;
            }
            if (!parse_double(argv[i + 1], muA_stats.mean_raw)) {
                std::cout << "ERROR INVALID --muA VALUE!" << std::endl;
                return 1;
            }
            muA_stats.provided = true;
            muA_stats.from_series = false;
            muA_stats.stderr_raw = 0.0;
            muA_stats.count = 1;
            muA_stats.source = "scalar";
            muA_scalar_seen = true;
            ++i;
        } else if (arg.compare(0, 6, "--muA=") == 0) {
            if (!parse_double(arg.substr(6), muA_stats.mean_raw)) {
                std::cout << "ERROR INVALID --muA VALUE!" << std::endl;
                return 1;
            }
            muA_stats.provided = true;
            muA_stats.from_series = false;
            muA_stats.stderr_raw = 0.0;
            muA_stats.count = 1;
            muA_stats.source = "scalar";
            muA_scalar_seen = true;
        } else if (arg == "--muA-series") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --muA-series REQUIRES A FILE!" << std::endl;
                return 1;
            }
            if (!load_muA_series(argv[i + 1], muA_stats)) {
                std::cout << "ERROR INVALID --muA-series FILE!" << std::endl;
                return 1;
            }
            muA_series_seen = true;
            ++i;
        } else if (arg.compare(0, 13, "--muA-series=") == 0) {
            if (!load_muA_series(arg.substr(13), muA_stats)) {
                std::cout << "ERROR INVALID --muA-series FILE!" << std::endl;
                return 1;
            }
            muA_series_seen = true;
        } else if (arg == "--muA-lammps-input") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --muA-lammps-input REQUIRES A FILE!" << std::endl;
                return 1;
            }
            muA_lammps_options.requested = true;
            muA_lammps_options.input_script = argv[i + 1];
            muA_lammps_seen = true;
            ++i;
        } else if (arg.compare(0, 19, "--muA-lammps-input=") == 0) {
            muA_lammps_options.requested = true;
            muA_lammps_options.input_script = arg.substr(19);
            muA_lammps_seen = true;
        } else if (arg == "--muA-lammps-lib") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --muA-lammps-lib REQUIRES A PATH!" << std::endl;
                return 1;
            }
            muA_lammps_options.library_path = argv[i + 1];
            ++i;
        } else if (arg.compare(0, 17, "--muA-lammps-lib=") == 0) {
            muA_lammps_options.library_path = arg.substr(17);
        } else if (arg == "--muA-lammps-numdiff") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --muA-lammps-numdiff REQUIRES A DELTA!" << std::endl;
                return 1;
            }
            if (!parse_double(argv[i + 1], muA_lammps_options.numdiff_delta) || !(muA_lammps_options.numdiff_delta > 0.0)) {
                std::cout << "ERROR INVALID --muA-lammps-numdiff DELTA!" << std::endl;
                return 1;
            }
            muA_lammps_options.use_numdiff = true;
            ++i;
        } else if (arg.compare(0, 21, "--muA-lammps-numdiff=") == 0) {
            if (!parse_double(arg.substr(21), muA_lammps_options.numdiff_delta) || !(muA_lammps_options.numdiff_delta > 0.0)) {
                std::cout << "ERROR INVALID --muA-lammps-numdiff DELTA!" << std::endl;
                return 1;
            }
            muA_lammps_options.use_numdiff = true;
        } else if (arg == "--muA-lammps-channel") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --muA-lammps-channel REQUIRES A CHANNEL!" << std::endl;
                return 1;
            }
            muA_lammps_options.channel = argv[i + 1];
            ++i;
        } else if (arg.compare(0, 21, "--muA-lammps-channel=") == 0) {
            muA_lammps_options.channel = arg.substr(21);
        } else if (arg == "--transform") {
            if (i + 1 >= argc) {
                std::cout << "ERROR --transform REQUIRES A METHOD!" << std::endl;
                return 1;
            }
            transform_method = normalize_transform(argv[i + 1]);
            ++i;
        } else if (arg.compare(0, 12, "--transform=") == 0) {
            transform_method = normalize_transform(arg.substr(12));
        } else if (!have_gtcut) {
            if (!parse_int(arg, Gtcut) || Gtcut < 0) {
                std::cout << "ERROR GTCUT!" << std::endl;
                return 1;
            }
            have_gtcut = true;
        } else if (!have_tail) {
            if (!parse_double(arg, tail_param)) {
                std::cout << "ERROR TAIL_PARAM!" << std::endl;
                return 1;
            }
            have_tail = true;
        } else {
            std::cout << "ERROR UNRECOGNIZED EXTRA ARGUMENT: " << arg << std::endl;
            return 1;
        }
    }

    const int muA_source_count = (muA_scalar_seen ? 1 : 0) + (muA_series_seen ? 1 : 0) + (muA_lammps_seen ? 1 : 0);
    if (muA_source_count > 1) {
        std::cout << "ERROR USE ONLY ONE muA SOURCE: --muA OR --muA-series OR --muA-lammps-input!" << std::endl;
        return 1;
    }
    if ((muA_lammps_options.use_numdiff || !muA_lammps_options.library_path.empty() || muA_lammps_options.channel != "auto") && !muA_lammps_seen) {
        std::cout << "ERROR --muA-lammps-lib / --muA-lammps-numdiff / --muA-lammps-channel REQUIRE --muA-lammps-input!" << std::endl;
        return 1;
    }

    if (dt <= 0.0 || !std::isfinite(dt)) {
        std::cout << "ERROR DT!" << std::endl;
        return 1;
    }
    if (mode != 0 && mode != 1) {
        std::cout << "ERROR MODE SET!" << std::endl;
        return 1;
    }
    if (!(tail_param > 0.0) || !std::isfinite(tail_param)) {
        std::cout << "ERROR TAIL_PARAM!" << std::endl;
        return 1;
    }
    if (transform_method != "direct" && transform_method != "irheo" && transform_method != "schwarzl" && transform_method != "caft") {
        std::cout << "ERROR UNKNOWN TRANSFORM METHOD!" << std::endl;
        return 1;
    }

    if (muA_lammps_seen) {
        if (!compute_muA_from_lammps(muA_lammps_options, mode, muA_lammps_result)) {
            std::cout << "ERROR LIBLAMMPS muA BACKEND FAILED: " << muA_lammps_result.error << std::endl;
            return 1;
        }
        muA_stats.provided = true;
        muA_stats.from_series = false;
        muA_stats.mean_raw = muA_lammps_result.muA_raw;
        muA_stats.stderr_raw = 0.0;
        muA_stats.count = 1;
        muA_stats.source = std::string("lammps:") + muA_lammps_options.input_script;
    }

    std::ifstream fin(input_file);
    std::ofstream fout(output_gt_file);
    std::ofstream fout2(output_gw_file);
    if (!fin.is_open() || fin.bad()) {
        std::cout << "ERROR INPUT FILE!" << std::endl;
        return 1;
    }
    if (!fout.is_open() || fout.bad() || !fout2.is_open() || fout2.bad()) {
        std::cout << "ERROR OUTPUT FILE!" << std::endl;
        return 1;
    }

    Correlator c;
    c.setsize(32, 16, 2);

    double valab, valac, valbc, valaa, valbb, valcc;
    int tt = 0;
    long line_number = 0;
    std::vector<double> Gt_reduced;
    std::vector<double> times_physical;
    const char *tensor_mode = (mode == 0) ? "shear3" : "tensor6";

    if (mode == 0) {
        c.initialize_mode0();
        while (fin >> valab >> valac >> valbc) {
            ++line_number;
            if (!std::isfinite(valab) || !std::isfinite(valac) || !std::isfinite(valbc)) {
                std::cout << "ERROR NON-FINITE INPUT AT LINE " << line_number << "!" << std::endl;
                return 1;
            }
            c.add_mode0(valab, valac, valbc);
            ++tt;
        }
        if (!fin.eof()) {
            std::cout << "ERROR INVALID INPUT FORMAT AT LINE " << (line_number + 1) << "!" << std::endl;
            return 1;
        }
        if (tt == 0) {
            std::cout << "ERROR EMPTY INPUT!" << std::endl;
            return 1;
        }
        c.evaluate_mode0(true);
        Gt_reduced.resize(c.npcorr, 0.0);
        for (unsigned int i = 0; i < c.npcorr; ++i) {
            Gt_reduced[i] = (c.f1_mode0[i] + c.f2_mode0[i] + c.f3_mode0[i]) / 3.0;
        }
    } else {
        c.initialize_mode1();
        while (fin >> valaa >> valbb >> valcc >> valab >> valac >> valbc) {
            ++line_number;
            if (!std::isfinite(valaa) || !std::isfinite(valbb) || !std::isfinite(valcc) ||
                !std::isfinite(valab) || !std::isfinite(valac) || !std::isfinite(valbc)) {
                std::cout << "ERROR NON-FINITE INPUT AT LINE " << line_number << "!" << std::endl;
                return 1;
            }
            c.add_mode1(valaa - valbb, valaa - valcc, valbb - valcc, valab, valac, valbc);
            ++tt;
        }
        if (!fin.eof()) {
            std::cout << "ERROR INVALID INPUT FORMAT AT LINE " << (line_number + 1) << "!" << std::endl;
            return 1;
        }
        if (tt == 0) {
            std::cout << "ERROR EMPTY INPUT!" << std::endl;
            return 1;
        }
        c.evaluate_mode1(true);
        Gt_reduced.resize(c.npcorr, 0.0);
        for (unsigned int i = 0; i < c.npcorr; ++i) {
            Gt_reduced[i] = (c.f1_mode1[i] + c.f2_mode1[i] + c.f3_mode1[i]) / 30.0
                           + (c.f4_mode1[i] + c.f5_mode1[i] + c.f6_mode1[i]) / 5.0;
        }
    }
    fin.close();

    if (c.npcorr == 0 || Gt_reduced.empty()) {
        std::cout << "ERROR NO CORRELATION POINTS GENERATED!" << std::endl;
        return 1;
    }

    times_physical.resize(c.npcorr, 0.0);
    for (unsigned int i = 0; i < c.npcorr; ++i) {
        times_physical[i] = c.t[i] * dt;
    }

    const char *branch = muA_stats.provided ? "absolute" : "reduced";
    const double muF_raw = Gt_reduced[0];
    const double Ge_raw = muA_stats.provided ? (muA_stats.mean_raw - muF_raw) : 0.0;
    const int effective_Gtcut = (Gtcut == -1) ? static_cast<int>(c.npcorr) : Gtcut;

    fout << "#dt=" << dt
         << " mode=" << mode
         << " tensor_mode=" << tensor_mode
         << " branch=" << branch
         << " centered=1"
         << " transform=" << transform_method
         << " Gtcut=" << effective_Gtcut
         << " tail_param=" << tail_param
         << std::endl;
    fout << "#raw_output=1 note=For_LJ_or_reduced_units_use_V_over_kBT;_for_LAMMPS_pressure-unit_input_use_V_over_(kBT*nktv2p)" << std::endl;
    fout << "#muF_raw=" << muF_raw;
    if (muA_stats.provided) {
        fout << " muA_raw=" << muA_stats.mean_raw
             << " muA_stderr_raw=" << muA_stats.stderr_raw
             << " muA_n=" << muA_stats.count
             << " muA_source=" << muA_stats.source
             << " Ge_raw=" << Ge_raw;
        if (muA_lammps_seen) {
            fout << " muA_lammps_units=" << muA_lammps_result.units
                 << " muA_lammps_channel=" << muA_lammps_result.channel
                 << " muA_lammps_library=" << muA_lammps_result.library_path
                 << " muA_phys=" << muA_lammps_result.muA_phys
                 << " muA_born_phys=" << muA_lammps_result.born_phys
                 << " muA_kinetic_phys=" << muA_lammps_result.kinetic_phys
                 << " muA_boltz=" << muA_lammps_result.boltz
                 << " muA_nktv2p=" << muA_lammps_result.nktv2p
                 << " muA_temperature=" << muA_lammps_result.temperature
                 << " muA_volume=" << muA_lammps_result.volume
                 << " muA_natoms=" << muA_lammps_result.natoms
                 << " muA_phys_to_raw_factor=" << muA_lammps_result.phys_to_raw_factor
                 << " muA_raw_to_phys_factor=" << muA_lammps_result.raw_to_phys_factor
                 << " muA_numdiff=" << (muA_lammps_result.used_numdiff ? 1 : 0);
        }
    }
    fout << std::endl;
    fout << "#Time " << (muA_stats.provided ? "G_abs_raw" : "G_red_raw") << std::endl;

    for (unsigned int i = 0; i < c.npcorr; ++i) {
        c.Gt[i] = Gt_reduced[i];
        const double output_gt = muA_stats.provided ? (Gt_reduced[i] + Ge_raw) : Gt_reduced[i];
        fout << times_physical[i] << " " << output_gt << std::endl;
    }
    fout.close();

    TrustedWindowInfo trusted_info = detect_trusted_window(times_physical, Gt_reduced, effective_Gtcut, tail_param);
    TransformResult tr;

    if (transform_method == "direct") {
        c.initialize_Gw();
        c.calculate_Gw(dt, tt, effective_Gtcut, tail_param);
        tr = collect_direct_result(c, dt, trusted_info);
    } else if (transform_method == "irheo") {
        tr = transform_irheo(times_physical, Gt_reduced, effective_Gtcut);
    } else if (transform_method == "schwarzl") {
        tr = transform_schwarzl(times_physical, Gt_reduced, effective_Gtcut);
    } else {
        tr = transform_caft_bounds(times_physical, Gt_reduced, dt, tt, effective_Gtcut, tail_param);
    }

    if (muA_stats.provided) {
        apply_absolute_shift(tr, Ge_raw);
    }

    fout2 << "#dt=" << dt
          << " mode=" << mode
          << " tensor_mode=" << tensor_mode
          << " branch=" << branch
          << " centered=1"
          << " transform=" << tr.method
          << " Gtcut=" << effective_Gtcut
          << " tail_param=" << tail_param
          << std::endl;
    fout2 << "#raw_output=1 trusted_time=" << tr.trusted_time
          << " trusted_points=" << tr.trusted_points
          << " tail_sigma=" << tr.tail_sigma
          << std::endl;
    fout2 << "#scaling_note=For_LJ_or_reduced_units_use_V_over_kBT;_for_LAMMPS_pressure-unit_input_use_V_over_(kBT*nktv2p)" << std::endl;
    fout2 << "#note=" << tr.note << std::endl;
    fout2 << "#muF_raw=" << muF_raw;
    if (muA_stats.provided) {
        fout2 << " muA_raw=" << muA_stats.mean_raw
              << " muA_stderr_raw=" << muA_stats.stderr_raw
              << " muA_n=" << muA_stats.count
              << " muA_source=" << muA_stats.source
              << " Ge_raw=" << Ge_raw;
        if (muA_lammps_seen) {
            fout2 << " muA_lammps_units=" << muA_lammps_result.units
                  << " muA_lammps_channel=" << muA_lammps_result.channel
                  << " muA_lammps_library=" << muA_lammps_result.library_path
                  << " muA_phys=" << muA_lammps_result.muA_phys
                  << " muA_born_phys=" << muA_lammps_result.born_phys
                  << " muA_kinetic_phys=" << muA_lammps_result.kinetic_phys
                  << " muA_boltz=" << muA_lammps_result.boltz
                  << " muA_nktv2p=" << muA_lammps_result.nktv2p
                  << " muA_temperature=" << muA_lammps_result.temperature
                  << " muA_volume=" << muA_lammps_result.volume
                  << " muA_natoms=" << muA_lammps_result.natoms
                  << " muA_phys_to_raw_factor=" << muA_lammps_result.phys_to_raw_factor
                  << " muA_raw_to_phys_factor=" << muA_lammps_result.raw_to_phys_factor
                  << " muA_numdiff=" << (muA_lammps_result.used_numdiff ? 1 : 0);
        }
    }
    fout2 << std::endl;
    fout2 << "#Frequency(1/time_unit) "
          << (muA_stats.provided ? "Gw_storage_abs_raw " : "Gw_storage_red_raw ")
          << "Gw_loss_raw Gw_mag_raw Gw_storage_low_raw Gw_storage_high_raw Gw_loss_low_raw Gw_loss_high_raw" << std::endl;

    for (std::size_t i = 0; i < tr.omega.size(); ++i) {
        fout2 << tr.omega[i] << " "
              << tr.storage[i] << " "
              << tr.loss[i] << " "
              << tr.magnitude[i] << " "
              << tr.storage_low[i] << " "
              << tr.storage_high[i] << " "
              << tr.loss_low[i] << " "
              << tr.loss_high[i] << std::endl;
    }
    fout2.close();

    return 0;
}
