#ifndef VELA_LAMMPS_MUA_BACKEND_H
#define VELA_LAMMPS_MUA_BACKEND_H

#include <string>

struct LammpsMuAOptions {
    bool requested;
    std::string input_script;
    std::string library_path;
    bool use_numdiff;
    double numdiff_delta;
    std::string channel;

    LammpsMuAOptions()
        : requested(false), input_script(), library_path(), use_numdiff(false), numdiff_delta(1.0e-5), channel("auto") {}
};

struct LammpsMuAResult {
    bool success;
    double muA_raw;
    double muA_phys;
    double born_phys;
    double kinetic_phys;
    double natoms;
    double temperature;
    double volume;
    double boltz;
    double nktv2p;
    double phys_to_raw_factor;
    double raw_to_phys_factor;
    std::string units;
    std::string channel;
    std::string library_path;
    bool used_numdiff;
    std::string note;
    std::string error;

    LammpsMuAResult()
        : success(false), muA_raw(0.0), muA_phys(0.0), born_phys(0.0), kinetic_phys(0.0),
          natoms(0.0), temperature(0.0), volume(0.0), boltz(0.0), nktv2p(1.0), phys_to_raw_factor(0.0), raw_to_phys_factor(0.0), units("unknown"),
          channel("auto"), library_path(""), used_numdiff(false), note(""), error("") {}
};

bool compute_muA_from_lammps(const LammpsMuAOptions& options, int mode, LammpsMuAResult& result);

#endif
