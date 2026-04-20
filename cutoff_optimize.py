#!/usr/bin/env python3
import argparse
import math
import subprocess
from pathlib import Path
from datetime import datetime

def read_qct(path: Path):
    data = {}
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            jf, cs = line.split()[:2]
            data[int(float(jf))] = float(cs)
    return data

def read_quantum(path: Path):
    data = {}
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            jf, cs = line.split()[:2]
            data[int(float(jf))] = float(cs)
    return data

def calc_metrics(qct, q):
    common = sorted(set(qct) & set(q))
    N = len(common)

    sum_qct = sum(qct[j] for j in common)
    sum_q   = sum(q[j]   for j in common)
    dInt = (sum_qct - sum_q) / sum_q

    sse = 0.0
    sae = 0.0
    for j in common:
        d = qct[j] - q[j]
        sse += d*d
        sae += abs(d)

    rmse_abs = math.sqrt(sse / N)
    l1_abs   = sae / N

    qmax = max(q[j] for j in common)
    eps  = 1e-3 * qmax
    sse_rel = sum(((qct[j] - q[j]) / (q[j] + eps))**2 for j in common)
    rmse_rel = math.sqrt(sse_rel / N)

    jf_peak_q    = max(common, key=lambda j: q[j])
    jf_peak_qct  = max(common, key=lambda j: qct[j])

    return {
        "dInt": dInt,
        "RMSE_abs": rmse_abs,
        "RMSE_rel": rmse_rel,
        "L1_abs": l1_abs,
        "sum_qct": sum_qct,
        "sum_q": sum_q,
        "peak_q": q[jf_peak_q],
        "peak_qct": qct[jf_peak_qct],
        "peak_ratio": qct[jf_peak_qct] / q[jf_peak_q],
        "jf_peak_q": jf_peak_q,
        "jf_peak_qct": jf_peak_qct,
        "N_common": N
    }

def linspace(a, b, n):
    if n == 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]

def frange(a, b, step):
    vals = []
    x = a
    while x <= b + 0.5 * step:
        vals.append(x)
        x += step
    return vals

def main():
    ap = argparse.ArgumentParser("Sweep Gaussian truncation cutoff")

    ap.add_argument("--binnew", default="./binnew_new")
    ap.add_argument("--traj", required=True)
    ap.add_argument("--dv", required=True, help="Δv (second argument to binnew_new)")
    ap.add_argument("--method", default="g")
    ap.add_argument("--vwidth", required=True)
    ap.add_argument("--jwidth", required=True)
    ap.add_argument("--quantum", required=True)

    ap.add_argument("--sigma-min", type=float, default=3.0)
    ap.add_argument("--sigma-max", type=float, default=3.6)
    ap.add_argument("--npoints", type=int, default=100)
    ap.add_argument("--step", type=float)
    ap.add_argument("--decimals", type=int, default=3)
    ap.add_argument("--outdir")
    ap.add_argument("--energy", required=True, help="Energy label for output filenames, e.g. 577")
    ap.add_argument("--initial-state", required=True, help="Initial state label for output filenames, e.g. v1j0")

    args = ap.parse_args()

    outdir = Path(args.outdir or f"cutoff_sweep_{datetime.now():%Y%m%d_%H%M%S}")
    outdir.mkdir(parents=True, exist_ok=True)

    quantum = read_quantum(Path(args.quantum))

    cutoffs = (
        frange(args.sigma_min, args.sigma_max, args.step)
        if args.step else
        linspace(args.sigma_min, args.sigma_max, args.npoints)
    )

    summary = outdir / "summary.tsv"
    with summary.open("w") as f:
        f.write(
            "# cutoff_sigma\tdInt\tRMSE_abs\tRMSE_rel\tL1_abs\t"
            "sum_qct\tsum_q\tpeak_qct\tpeak_q\tpeak_ratio\t"
            "jf_peak_qct\tjf_peak_q\tN_common\tqct_file\n"
        )

    for i, n in enumerate(cutoffs, 1):
        cutoff = f"{n:.{args.decimals}f}"
        qct_out = outdir / f"e{args.energy}_v{args.initial_state}_dv{args.dv}_{cutoff}sig.cbt"

        cmd = [
            args.binnew,
            args.traj,
            args.dv,
            args.method,
            args.vwidth,
            args.jwidth,
            cutoff
        ]

        print(f"[{i}/{len(cutoffs)}] σ={cutoff}")
        with qct_out.open("w") as f:
            subprocess.run(cmd, stdout=f, check=True)

        qct = read_qct(qct_out)
        m = calc_metrics(qct, quantum)

        with summary.open("a") as f:
            f.write(
                f"{cutoff}\t{m['dInt']:.6g}\t{m['RMSE_abs']:.6g}\t"
                f"{m['RMSE_rel']:.6g}\t{m['L1_abs']:.6g}\t"
                f"{m['sum_qct']:.6g}\t{m['sum_q']:.6g}\t"
                f"{m['peak_qct']:.6g}\t{m['peak_q']:.6g}\t"
                f"{m['peak_ratio']:.6g}\t"
                f"{m['jf_peak_qct']}\t{m['jf_peak_q']}\t"
                f"{m['N_common']}\t{qct_out.name}\n"
            )

    print(f"\nDone. Summary → {summary}")

if __name__ == "__main__":
    main()

