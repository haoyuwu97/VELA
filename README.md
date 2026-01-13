# Modulus-time and modulus-frequency relationships (G(t) and G($\omega$)) based on autocorrelation function

## Overview
This package computes stress autocorrelation functions in the time domain and transforms them to modulus-frequency relationships. The raw output should be scaled by `V/kBT` to obtain physical moduli, as noted in the CLI help.

## Units and physical interpretation
The shear modulus from equilibrium stress autocorrelation is typically defined as:

```
G(t) = (V / kBT) * <σ(0) σ(t)>
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

## Output format
- **G(t)**: `time  Gt`
- **G(ω)**: `frequency  Gw_storage  Gw_loss  Gw`

Output headers include metadata for reproducibility:

```
#dt=... mode=... Gtcut=... tail_param=...
```

## Validation example (Ornstein-Uhlenbeck)
A small OU-process generator is provided to validate that the computed `G(t)` follows an exponential decay.

1. Generate synthetic mode-0 data:
   ```bash
   python3 scripts/generate_ou.py --nsteps 5000 --dt 0.01 --tau 1.0 --sigma 1.0
   ```
2. Build and run:
   ```bash
   make
   ./Gt_Gw data/ou_mode0.dat out_gt.dat out_gw.dat 0.01 0
   ```
3. Compare `out_gt.dat` to `data/expected_gt.csv` (they should show the same exponential trend).
