#include "lammps_mua_backend.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>

#if defined(__linux__) || defined(__APPLE__)
#include <dlfcn.h>
#include <unistd.h>
#endif

namespace {

typedef void* (*LammpsOpenNoMpiFn)(int, char **, void *);
typedef void (*LammpsCloseFn)(void *);
typedef void (*LammpsFileFn)(void *, const char *);
typedef char* (*LammpsCommandFn)(void *, const char *);
typedef double (*LammpsGetNatomsFn)(void *);
typedef double (*LammpsGetThermoFn)(void *, const char *);
typedef void* (*LammpsExtractGlobalFn)(void *, const char *);
typedef int (*LammpsHasErrorFn)(void *);
typedef int (*LammpsGetLastErrorMessageFn)(void *, char *, int);
typedef int (*LammpsSetShowErrorFn)(void *, int);

struct DynamicLib {
    void *handle;
    DynamicLib() : handle(NULL) {}
    ~DynamicLib() {
#if defined(__linux__) || defined(__APPLE__)
        if (handle) dlclose(handle);
#endif
    }
};

struct LammpsApi {
    DynamicLib lib;
    LammpsOpenNoMpiFn open_no_mpi;
    LammpsCloseFn close;
    LammpsFileFn file;
    LammpsCommandFn command;
    LammpsGetNatomsFn get_natoms;
    LammpsGetThermoFn get_thermo;
    LammpsExtractGlobalFn extract_global;
    LammpsHasErrorFn has_error;
    LammpsGetLastErrorMessageFn get_last_error_message;
    LammpsSetShowErrorFn set_show_error;
    std::string loaded_path;

    LammpsApi()
        : open_no_mpi(NULL), close(NULL), file(NULL), command(NULL), get_natoms(NULL), get_thermo(NULL),
          extract_global(NULL), has_error(NULL), get_last_error_message(NULL), set_show_error(NULL), loaded_path() {}
};

bool parse_numbers_line(const std::string &line, std::vector<double> &vals) {
    vals.clear();
    std::stringstream ss(line);
    double x = 0.0;
    while (ss >> x) vals.push_back(x);
    return !vals.empty();
}

std::string lowercase_copy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string make_tmpfile_path() {
#if defined(__linux__) || defined(__APPLE__)
    char templ[] = "/tmp/vela_bornXXXXXX";
    int fd = mkstemp(templ);
    if (fd >= 0) {
        close(fd);
        std::remove(templ);
        return std::string(templ);
    }
#endif
    return std::string("vela_born_tmp.out");
}

bool file_exists(const std::string &path) {
    std::ifstream fin(path.c_str());
    return fin.good();
}

void *load_symbol(void *lib, const char *name) {
#if defined(__linux__) || defined(__APPLE__)
    return dlsym(lib, name);
#else
    (void)lib;
    (void)name;
    return NULL;
#endif
}

bool load_lammps_api(const LammpsMuAOptions &options, LammpsApi &api, std::string &err) {
#if !defined(__linux__) && !defined(__APPLE__)
    err = "dynamic liblammps loading is currently implemented only for POSIX platforms in this VELA build";
    return false;
#else
    std::vector<std::string> candidates;
    if (!options.library_path.empty()) {
        candidates.push_back(options.library_path);
    } else {
        candidates.push_back("liblammps.so");
        candidates.push_back("liblammps_mpi.so");
        candidates.push_back("liblammps.dylib");
        candidates.push_back("liblammps_mpi.dylib");
        candidates.push_back("liblammps.dll");
        candidates.push_back("lammps.dll");
    }
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        api.lib.handle = dlopen(candidates[i].c_str(), RTLD_LAZY | RTLD_LOCAL);
        if (api.lib.handle) {
            api.loaded_path = candidates[i];
            break;
        }
    }
    if (!api.lib.handle) {
        err = "failed to load a liblammps shared library; pass --muA-lammps-lib PATH or make liblammps discoverable";
        return false;
    }

    api.open_no_mpi = reinterpret_cast<LammpsOpenNoMpiFn>(load_symbol(api.lib.handle, "lammps_open_no_mpi"));
    api.close = reinterpret_cast<LammpsCloseFn>(load_symbol(api.lib.handle, "lammps_close"));
    api.file = reinterpret_cast<LammpsFileFn>(load_symbol(api.lib.handle, "lammps_file"));
    api.command = reinterpret_cast<LammpsCommandFn>(load_symbol(api.lib.handle, "lammps_command"));
    api.get_natoms = reinterpret_cast<LammpsGetNatomsFn>(load_symbol(api.lib.handle, "lammps_get_natoms"));
    api.get_thermo = reinterpret_cast<LammpsGetThermoFn>(load_symbol(api.lib.handle, "lammps_get_thermo"));
    api.extract_global = reinterpret_cast<LammpsExtractGlobalFn>(load_symbol(api.lib.handle, "lammps_extract_global"));
    api.has_error = reinterpret_cast<LammpsHasErrorFn>(load_symbol(api.lib.handle, "lammps_has_error"));
    api.get_last_error_message = reinterpret_cast<LammpsGetLastErrorMessageFn>(load_symbol(api.lib.handle, "lammps_get_last_error_message"));
    api.set_show_error = reinterpret_cast<LammpsSetShowErrorFn>(load_symbol(api.lib.handle, "lammps_set_show_error"));

    if (!api.open_no_mpi || !api.close || !api.file || !api.command || !api.get_natoms || !api.get_thermo || !api.extract_global) {
        err = "loaded liblammps is missing one or more required C-library symbols";
        return false;
    }
    return true;
#endif
}

std::string fetch_last_error(const LammpsApi &api, void *handle) {
    if (!api.has_error || !api.get_last_error_message) return std::string();
    if (!api.has_error(handle)) return std::string();
    char buffer[1024];
    std::memset(buffer, 0, sizeof(buffer));
    api.get_last_error_message(handle, buffer, static_cast<int>(sizeof(buffer)));
    return std::string(buffer);
}

bool run_lammps_file(const LammpsApi &api, void *handle, const std::string &path, std::string &err) {
    api.file(handle, path.c_str());
    err = fetch_last_error(api, handle);
    return err.empty();
}

bool run_lammps_command(const LammpsApi &api, void *handle, const std::string &cmd, std::string &err) {
    api.command(handle, cmd.c_str());
    err = fetch_last_error(api, handle);
    return err.empty();
}

double get_global_double(const LammpsApi &api, void *handle, const char *name, bool &ok) {
    ok = false;
    void *ptr = api.extract_global(handle, name);
    if (!ptr) return 0.0;
    ok = true;
    return *reinterpret_cast<double *>(ptr);
}

std::string get_global_string(const LammpsApi &api, void *handle, const char *name) {
    void *ptr = api.extract_global(handle, name);
    if (!ptr) return std::string();
    return std::string(reinterpret_cast<const char *>(ptr));
}

double compute_mode0_muA_phys(const std::vector<double> &born, double vol, double nktv2p, double kin_phys, const std::string &channel) {
    const double c44 = born[3] / vol * nktv2p;
    const double c55 = born[4] / vol * nktv2p;
    const double c66 = born[5] / vol * nktv2p;
    if (channel == "yz") return c44 + kin_phys;
    if (channel == "xz") return c55 + kin_phys;
    if (channel == "xy") return c66 + kin_phys;
    return (c44 + c55 + c66) / 3.0 + kin_phys;
}

double compute_mode1_muA_phys(const std::vector<double> &born, double vol, double nktv2p, double kin_phys, const std::string &channel) {
    const double c11 = born[0] / vol * nktv2p;
    const double c22 = born[1] / vol * nktv2p;
    const double c33 = born[2] / vol * nktv2p;
    const double c44 = born[3] / vol * nktv2p;
    const double c55 = born[4] / vol * nktv2p;
    const double c66 = born[5] / vol * nktv2p;
    const double c12 = born[6] / vol * nktv2p;
    const double c13 = born[7] / vol * nktv2p;
    const double c23 = born[11] / vol * nktv2p;

    if (channel == "yz") return c44 + kin_phys;
    if (channel == "xz") return c55 + kin_phys;
    if (channel == "xy") return c66 + kin_phys;
    if (channel == "avg3") return (c44 + c55 + c66) / 3.0 + kin_phys;

    const double born_iso = ((c11 + c22 + c33) - (c12 + c13 + c23)) / 15.0 + (c44 + c55 + c66) / 5.0;
    const double kin_iso = 0.8 * kin_phys;
    return born_iso + kin_iso;
}

} // namespace

bool compute_muA_from_lammps(const LammpsMuAOptions& options, int mode, LammpsMuAResult& result) {
    result = LammpsMuAResult();
    result.channel = options.channel;
    result.used_numdiff = options.use_numdiff;

    if (!options.requested) {
        result.error = "lammps muA backend requested=false";
        return false;
    }
    if (options.input_script.empty()) {
        result.error = "--muA-lammps-input is required for the liblammps backend";
        return false;
    }
    if (!file_exists(options.input_script)) {
        result.error = "the specified --muA-lammps-input file does not exist";
        return false;
    }

    LammpsApi api;
    std::string err;
    if (!load_lammps_api(options, api, err)) {
        result.error = err;
        return false;
    }
    result.library_path = api.loaded_path;

    const char *args_raw[] = {"liblammps", "-nocite", "-log", "none", "-echo", "none", NULL};
    char **argv = const_cast<char **>(args_raw);
    const int argc = 6;
    void *handle = api.open_no_mpi(argc, argv, NULL);
    if (!handle) {
        result.error = "failed to create a LAMMPS instance via lammps_open_no_mpi";
        return false;
    }

    if (api.set_show_error) api.set_show_error(handle, 0);

    bool ok = true;
    do {
        if (!run_lammps_file(api, handle, options.input_script, err)) { ok = false; break; }
        if (!run_lammps_command(api, handle, "compute __vela_pvir all pressure NULL virial", err)) { ok = false; break; }
        std::stringstream ccmd;
        if (options.use_numdiff) {
            ccmd << "compute __vela_born all born/matrix numdiff " << options.numdiff_delta << " __vela_pvir";
        } else {
            ccmd << "compute __vela_born all born/matrix";
        }
        if (!run_lammps_command(api, handle, ccmd.str(), err)) { ok = false; break; }
        if (!run_lammps_command(api, handle, "run 0 post no", err)) { ok = false; break; }

        const std::string tmpfile = make_tmpfile_path();
        for (int i = 1; i <= 21; ++i) {
            std::stringstream vcmd;
            vcmd << "variable __vela_b" << i << " equal c___vela_born[" << i << "]";
            if (!run_lammps_command(api, handle, vcmd.str(), err)) { ok = false; break; }
        }
        if (!ok) break;
        std::stringstream pcmd;
        pcmd << "print \"";
        for (int i = 1; i <= 21; ++i) {
            if (i > 1) pcmd << ' ';
            pcmd << "${__vela_b" << i << "}";
        }
        pcmd << "\" file " << tmpfile << " screen no";
        if (!run_lammps_command(api, handle, pcmd.str(), err)) { ok = false; break; }

        std::ifstream fin(tmpfile.c_str());
        std::string line;
        if (!std::getline(fin, line)) {
            err = "failed to read born matrix values from the temporary liblammps output file";
            ok = false;
            std::remove(tmpfile.c_str());
            break;
        }
        std::remove(tmpfile.c_str());

        std::vector<double> born;
        if (!parse_numbers_line(line, born) || born.size() < 21) {
            err = "invalid born matrix output received from liblammps";
            ok = false;
            break;
        }
        born.resize(21);

        result.natoms = api.get_natoms(handle);
        result.temperature = api.get_thermo(handle, "temp");
        result.volume = api.get_thermo(handle, "vol");
        result.units = get_global_string(api, handle, "units");
        const std::string units_lower = lowercase_copy(result.units);

        bool got_boltz = false;
        bool got_nktv2p = false;
        result.boltz = get_global_double(api, handle, "boltz", got_boltz);
        result.nktv2p = get_global_double(api, handle, "nktv2p", got_nktv2p);
        if (!got_boltz) result.boltz = 1.0;
        if (!got_nktv2p) result.nktv2p = 1.0;
        if (!(result.volume > 0.0) || !(result.boltz > 0.0) || !(result.temperature > 0.0)) {
            err = "liblammps returned invalid thermodynamic state (temperature, volume, or boltz <= 0) for conversion to VELA raw units";
            ok = false;
            break;
        }
        if (units_lower != "lj" && !units_lower.empty() && !got_nktv2p) {
            err = "liblammps did not expose nktv2p for a non-lj unit style; raw/physical conversion would be ambiguous";
            ok = false;
            break;
        }

        const double pressure_scale = result.boltz * result.temperature / result.volume * result.nktv2p;
        if (!(pressure_scale > 0.0) || !std::isfinite(pressure_scale)) {
            err = "failed to construct a valid thermal pressure scale kB*T*nktv2p/V for VELA raw-unit conversion";
            ok = false;
            break;
        }
        result.phys_to_raw_factor = pressure_scale;
        result.raw_to_phys_factor = 1.0 / pressure_scale;

        const double kin_phys = result.natoms * pressure_scale;
        result.kinetic_phys = kin_phys;

        std::string channel = options.channel;
        for (std::size_t i = 0; i < channel.size(); ++i) channel[i] = static_cast<char>(std::tolower(channel[i]));
        if (channel.empty() || channel == "auto") {
            channel = (mode == 0) ? "avg3" : "iso6";
        }
        if (!(channel == "xy" || channel == "xz" || channel == "yz" || channel == "avg3" || channel == "iso6")) {
            err = "unsupported --muA-lammps-channel; use auto|xy|xz|yz|avg3|iso6";
            ok = false;
            break;
        }
        if (mode == 0 && channel == "iso6") channel = "avg3";
        result.channel = channel;

        if (mode == 0) {
            result.muA_phys = compute_mode0_muA_phys(born, result.volume, result.nktv2p, kin_phys, channel);
            result.born_phys = result.muA_phys - kin_phys;
        } else {
            result.muA_phys = compute_mode1_muA_phys(born, result.volume, result.nktv2p, kin_phys, channel);
            if (channel == "yz" || channel == "xz" || channel == "xy" || channel == "avg3") {
                result.born_phys = result.muA_phys - kin_phys;
            } else {
                result.born_phys = result.muA_phys - 0.8 * kin_phys;
            }
        }

        result.muA_raw = result.muA_phys * result.phys_to_raw_factor;
        result.note = "muA obtained from liblammps using compute born/matrix plus the kinetic affine term and then converted to VELA raw units via muA_raw = muA_phys * (kB*T*nktv2p / V); for lj-style reduced units nktv2p=1";
        result.success = true;
    } while (false);

    if (!ok) {
        result.error = err.empty() ? std::string("liblammps muA backend failed") : err;
    }

    api.close(handle);
    return result.success;
}
