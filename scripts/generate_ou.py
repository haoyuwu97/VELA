#!/usr/bin/env python3
import argparse
import math
import random
from pathlib import Path


def generate_ou(nsteps, dt, tau, sigma, seed):
    rng = random.Random(seed)
    alpha = math.exp(-dt / tau)
    noise_scale = sigma * math.sqrt(1.0 - alpha * alpha)
    x = 0.0
    series = []
    for _ in range(nsteps):
        x = alpha * x + noise_scale * rng.gauss(0.0, 1.0)
        series.append(x)
    return series


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nsteps", type=int, default=5000)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--tau", type=float, default=1.0)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--out", type=Path, default=Path("data/ou_mode0.dat"))
    parser.add_argument("--expected", type=Path, default=Path("data/expected_gt.csv"))
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.expected.parent.mkdir(parents=True, exist_ok=True)

    ab = generate_ou(args.nsteps, args.dt, args.tau, args.sigma, args.seed)
    ac = generate_ou(args.nsteps, args.dt, args.tau, args.sigma, args.seed + 1)
    bc = generate_ou(args.nsteps, args.dt, args.tau, args.sigma, args.seed + 2)

    with args.out.open("w", encoding="utf-8") as handle:
        for a, b, c in zip(ab, ac, bc):
            handle.write(f"{a} {b} {c}\n")

    with args.expected.open("w", encoding="utf-8") as handle:
        handle.write("time,expected_gt\n")
        for i in range(args.nsteps):
            t = i * args.dt
            expected = args.sigma * args.sigma * math.exp(-t / args.tau)
            handle.write(f"{t},{expected}\n")


if __name__ == "__main__":
    main()
