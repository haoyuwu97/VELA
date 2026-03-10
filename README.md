# VELA: ViscoElastic Linear-response Analyzer

## Overview

VELA computes equilibrium stress autocorrelation functions in the time domain and transforms them into modulus-frequency relationships from MD stress time series.

The current code distinguishes two physical branches automatically:

- **Reduced branch (default):** used automatically when only equilibrium stress components are supplied. The output is the centered fixed-strain stress autocorrelation, i.e. the reduced modulus quantity.
- **Absolute branch:** activated when a raw affine modulus is supplied via `--muA VALUE`, a statistically averaged `--muA-series FILE`, or the new `--muA-lammps-input FILE` backend.

> **Scaling convention:** VELA still outputs **raw** moduli. For LJ/CG-style reduced units (`nktv2p = 1`), multiply every reported modulus and `muA` by `V/(k_B T)` to obtain physical moduli. If the input stresses are in LAMMPS pressure units (`real`, `metal`, `si`, `cgs`, ...), the correct conversion is `V/(k_B T nktv2p)`, because the raw correlators carry pressure-squared units while LAMMPS pressure uses the extra conversion constant `nktv2p`.

---

## Physical Background

### Reduced branch

For stress-only fixed-strain equilibrium input, VELA reports the centered stress autocorrelation quantity

```math
G_{\mathrm{red}}(t) = \langle \delta\sigma(0)\,\delta\sigma(t)\rangle
```

in raw units. This is the physically safe default when no affine modulus information is supplied.

### Absolute branch

If a raw affine modulus `muA_raw` is supplied with `--muA` or `--muA-series`, VELA estimates the raw stress-fluctuation contribution from the zero-time value of the centered correlator,

```math
\mu_{F,\mathrm{raw}} = G_{\mathrm{red}}(0)
```

and then forms the raw plateau shift

```math
G_{e,\mathrm{raw}} = \mu_{A,\mathrm{raw}} - \mu_{F,\mathrm{raw}}.
```

The reported absolute raw modulus becomes

```math
G_{\mathrm{abs,raw}}(t)=G_{\mathrm{red,raw}}(t)+G_{e,\mathrm{raw}}.
```

In frequency space the same constant shift is added only to the storage modulus.

---

## Build

```bash
make
```

---

## Usage

```bash
./Gt_Gw INPUT_FILE OUTPUT_FILE_Gt OUTPUT_FILE_Gw DT MODE [GTCUT] [TAIL_PARAM] [--muA VALUE] [--muA-series FILE] [--muA-lammps-input FILE] [--muA-lammps-lib PATH] [--muA-lammps-numdiff DELTA] [--muA-lammps-channel CHANNEL] [--transform METHOD]
```

### Required arguments

- `INPUT_FILE`: stress components per time step.
- `OUTPUT_FILE_Gt`: time-domain output.
- `OUTPUT_FILE_Gw`: frequency-domain output.
- `DT`: time step of the smallest data spacing.
- `MODE`: `0` for `xy,xz,yz` (3 columns), `1` for `xx,yy,zz,xy,xz,yz` (6 columns).

### Optional arguments

- `GTCUT`: maximum time index used by transform-range selection or trusted-window logic (defaults to the full length).
- `TAIL_PARAM`: k-sigma threshold for the statistical confidence cutoff / trusted-window detector (default `2.0`).
- `--muA VALUE`: switches to the **absolute** branch using the supplied raw affine modulus.
- `--muA-series FILE`: switches to the **absolute** branch using the mean of a raw affine-modulus series. Each non-comment line may contain either one value or multiple columns; VELA uses the last numeric column on each line.
- `--muA-lammps-input FILE`: switches to the **absolute** branch by asking VELA to create a `liblammps` instance and compute `muA` from a complete LAMMPS setup script. The script must fully define the system, force field, and how the target configuration is read (`read_data`, `read_restart`, `read_dump`, etc.).
- `--muA-lammps-lib PATH`: optional path to the shared `liblammps` library. If omitted, VELA tries common library names such as `liblammps.so` and `liblammps_mpi.so`.
- `--muA-lammps-numdiff DELTA`: optional fallback to `compute born/matrix numdiff DELTA __vela_pvir` for cases where analytic Born derivatives are unavailable. The LAMMPS documentation recommends sensitivity testing and notes that typical useful deltas are often in the `1e-5` to `1e-6` range.
- `--muA-lammps-channel CHANNEL`: choose `auto`, `xy`, `xz`, `yz`, `avg3`, or `iso6`. `auto` means `avg3` for mode 0 and `iso6` for mode 1.
- `--transform METHOD`: choose `direct`, `irheo`, `schwarzl`, or `caft` (default: `direct`).

---

## Input Modes

### Mode 0 (3 components)

Each line of the input file has

```text
xy xz yz
```

VELA computes the centered average of the three shear autocorrelation functions.

### Mode 1 (6 components)

Each line of the input file has

```text
xx yy zz xy xz yz
```

VELA internally converts the normal stresses to normal-stress differences

```math
N_{xy}=\sigma_{xx}-\sigma_{yy},\qquad N_{xz}=\sigma_{xx}-\sigma_{zz},\qquad N_{yz}=\sigma_{yy}-\sigma_{zz}
```

and then evaluates the weighted isotropic average

```math
G_{\mathrm{red}}(t)=\frac{\langle N_{xy}(0)N_{xy}(t)\rangle+\langle N_{xz}(0)N_{xz}(t)\rangle+\langle N_{yz}(0)N_{yz}(t)\rangle}{30}
+\frac{\langle\sigma_{xy}(0)\sigma_{xy}(t)\rangle+\langle\sigma_{xz}(0)\sigma_{xz}(t)\rangle+\langle\sigma_{yz}(0)\sigma_{yz}(t)\rangle}{5}
```

with **centered** correlations used by default.

---

### `--muA-lammps-input FILE`

This backend is the current minimal autonomous `muA` route.

VELA dynamically loads the LAMMPS shared library and then:

1. reads a complete LAMMPS setup script with `lammps_file()`;
2. creates `compute __vela_born all born/matrix` (or the `numdiff` variant);
3. executes `run 0 post no`;
4. reads thermodynamic data (`temp`, `vol`, `natoms`) and LAMMPS unit constants (`boltz`, `nktv2p`);
5. forms an affine modulus in pressure units and then converts it to the current VELA raw convention through

```math
\mu_{A,\mathrm{raw}} = \mu_{A,\mathrm{phys}}\,\frac{k_B T\,nktv2p}{V}.
```

This reduces to the familiar LJ/CG expression when `nktv2p = 1`. When the `liblammps` backend is used, VELA also reports `muA_raw_to_phys_factor = V/(k_B T nktv2p)` in the output header so the same factor can be applied to every reported raw modulus.

For mode 0, the default channel is the three-shear average

```math
\mu_A^{(3\mathrm{shear})}=\frac{C_{44}+C_{55}+C_{66}}{3}.
```

For mode 1, the default channel is the isotropic scalar consistent with the current normal-difference/shear projection,

```math
\mu_A^{(\mathrm{iso6})}=\frac{(C_{11}+C_{22}+C_{33})-(C_{12}+C_{13}+C_{23})}{15}+\frac{C_{44}+C_{55}+C_{66}}{5}.
```

The kinetic affine term is included. For the six-component isotropic projection, this implies the corresponding projected kinetic contribution. The LAMMPS Born-matrix documentation and elastic-constants howto are the physical basis for this backend.

## Transform Methods

### `--transform direct`

The current VELA direct transform on the multi-tau grid. It preserves the old code path and interprets the old k-sigma tail cutoff as a trusted-window detector.

### `--transform irheo`

A baseline implementation of the non-oversampled i-RheoFT logic used in the open RepTate G(t) application.

### `--transform schwarzl`

A **Schwarzl-style numerical baseline** on a log-resampled time grid. The current implementation is intended as a transparent numerical reference and **not** as a byte-identical reproduction of RepTate's helper library.

### `--transform caft`

`CAFT-Bounds` (Constrained Analytic Fourier Transform with Bounds), VELA's current main method:

1. detect a trusted time window from the time-domain data;
2. fit a nonnegative multi-exponential tail continuation;
3. transform the trusted prefix plus the constrained tail;
4. report lower/upper envelopes for storage and loss moduli.

---

## Current Output Semantics

VELA writes explicit metadata in the output headers:

- `branch=reduced` when no `muA` information is provided.
- `branch=absolute` when `--muA`, `--muA-series`, or `--muA-lammps-input` is supplied.
- `centered=1` because centered correlations are now the scientific default.
- `transform=...` to identify the frequency-domain conversion route.

### Time-domain output

The first column is time. The second column is:

- reduced modulus in the reduced branch;
- absolute modulus in the absolute branch.

### Frequency-domain output

The first four columns are always

```text
frequency storage loss magnitude
```

and the next four columns are

```text
storage_low storage_high loss_low loss_high
```

For `direct`, `irheo`, and `schwarzl`, the bounds currently collapse to the point estimate.
For `caft`, the bounds reflect the trusted-window and constrained-tail inference.

---

## Files

- `main.cc`: CLI, automatic reduced/absolute branch selection, muA handling, output logic.
- `correlator.h`: shared declarations.
- `correlator_mode0.cc`: 3-component correlator.
- `correlator_mode1.cc`: 6-component correlator.
- `correlator_Gw.cc`: existing direct frequency-domain transform.
- `transform_methods.h/.cc`: i-RheoFT baseline, Schwarzl-style baseline, CAFT-Bounds, trusted-window logic.
- `lammps_mua_backend.h/.cc`: optional runtime `liblammps` backend for autonomous `muA` extraction.
- `Makefile`: build instructions.

---

## Final scope of this package version

- `--muA`, `--muA-series`, and `--muA-lammps-input` are all supported.
- The autonomous `muA` route is **by design** the `liblammps` route: the user must provide a **complete LAMMPS setup script** that already knows how to read the target `data` / `restart` / `dump` configuration and define the full force field. VELA then appends the Born/virial commands and executes `run 0`.
- There is **no separate native VELA force-field reimplementation** in this package version. This is intentional: the package closes around the `liblammps` backend rather than maintaining a duplicate force engine.
- `direct`, `irheo`, `schwarzl`, and `caft` are implemented.
- The present `schwarzl` transform is intentionally documented as a **Schwarzl-style numerical baseline**, not as an exact reproduction of RepTate's C helper.
- `CAFT-Bounds` is a **v1 implementation** using trusted-window detection plus nonnegative multi-exponential tail continuation.
- The `liblammps` backend requires a **strictly positive temperature** because the current VELA raw-unit convention is thermal and uses the prefactor `k_B T nktv2p / V`.
