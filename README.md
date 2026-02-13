# VELA: ViscoElastic Linear-response Analyzer

## Overview
This package computes stress autocorrelation functions in the time domain and transforms them to modulus-frequency relationships. The raw output should be scaled by `V/kBT` to obtain physical moduli, as noted in the CLI help.

## Units and physical interpretation
The shear modulus from equilibrium stress autocorrelation is typically defined as:
```math
G(t)=\frac{V}{kBT}\langle\sigma(0)\sigma(t)\rangle
```
where `σ` is a stress component, `V` is the system volume, `T` is temperature, and `kB` is Boltzmann's constant. Ensure that:
- Stress components in the input file match your MD engine's output units.
- `DT` matches the time unit of the input series.

**Example (LAMMPS):**
- If stresses are reported in pressure units (e.g., "real" -> atm, "metal" -> bar), then `σ` uses those units.
- Volume uses the same length units cubed as your simulation box.
- `kBT` uses the energy unit of the simulation.

Always verify unit consistency before scaling `G(t)` or `G(ω)`.

## Build
```bash
make
```

## Usage
```bash
./Gt_Gw INPUT_FILE OUTPUT_FILE_Gt OUTPUT_FILE_Gw DT MODE [GTCUT] [TAIL_PARAM]
```

**Arguments**
- `INPUT_FILE`: stress components per time step.
- `OUTPUT_FILE_Gt`: time-domain modulus output (two columns: time, G(t)).
- `OUTPUT_FILE_Gw`: frequency-domain modulus output.
- `DT`: time step of the smallest data spacing.
- `MODE`: `0` for `xy,xz,yz` (3 columns), `1` for `xx,yy,zz,xy,xz,yz` (6 columns).
- `GTCUT` (optional): maximum time index considered for tail handling (defaults to full length).
- `TAIL_PARAM` (optional): k-sigma threshold for the statistical confidence cutoff (default `2.0`).

## Input format
- **Mode 0**: each line has `xy xz yz`.
- **Mode 1**: each line has `xx yy zz xy xz yz` (normal components first).

**Mode 0 (3 components)**   
The time-domain modulus is computed from the average of the three shear autocorrelation functions:
```math
G(t)=\frac{V}{kBT}\frac{\langle\sigma_{xy}(0)\sigma_{xy}(t)\rangle+\langle\sigma_{xz}(0)\sigma_{xz}(t)\rangle+\langle\sigma_{yz}(0)\sigma_{yz}(t)\rangle}{3}
```
Use Mode 0 when your input file only contains shear stresses.

**Mode 1 (6 components)**  
The time-domain modulus is computed from both normal and shear components with weighting:
```math
G(t)=\frac{V}{kBT}\left(\frac{\langle N_{xy}(0)N_{xy}(t)\rangle+\langle N_{yz}(0)N_{yz}(t)\rangle+\langle N_{xz}(0)N_{xz}(t)\rangle}{30}+\frac{\langle\sigma_{xy}(0)\sigma_{xy}(t)\rangle+\langle\sigma_{yz}(0)\sigma_{yz}(t)\rangle+\langle\sigma_{xz}(0)\sigma_{xz}(t)\rangle}{5}\right)
```
```math
N_{\alpha\beta}(t)=\sigma_{\alpha\alpha}(t)-\sigma_{\beta\beta}(t)
```
Use Mode 1 when you have full stress tensor components and want the isotropic average.

## Output format
- **G(t)**: `time  Gt`
- **G(ω)**: `frequency  Gw_storage  Gw_loss  Gw`

Output headers include metadata for reproducibility:

```
#dt=... mode=... Gtcut=... tail_param=...
```
