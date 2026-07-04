#!/usr/bin/env python3
"""
Detailed-balance checker for rotationally inelastic dv = 0 data.

This version is designed for a *set* of input files such as

    v1j0ExpHP.rc  v1j2ExpHP.rc  v1j4ExpHP.rc ...

where each file contains rows

    j_f   sigma_if   err_if

For every available pair i <-> f it calculates the Excel-style detailed
balance comparison:

    forward_weight = (2*j_i + 1) * exp[-E_i/(k_B T)] * sigma(i -> f)
    reverse_weight = (2*j_f + 1) * exp[-E_f/(k_B T)] * sigma(f -> i)
    DB_ratio       = forward_weight / reverse_weight

If detailed balance is satisfied, DB_ratio should be close to 1.

Energy input file formats accepted:

    v1j0   379.98509
    v1j2   382.927109

or

    1   0   379.98509
    1   2   382.927109

Default energy unit is cm^-1, matching the spreadsheet approach.
"""

import argparse
import csv
import glob
import math
import os
import re
from collections import defaultdict
from statistics import mean, stdev

# Exact SI constants.  Energies in cm^-1 are converted to J by h*c*(cm/s).
KB = 1.380649e-23          # J K^-1
H = 6.62607015e-34         # J s
C_CM_S = 2.99792458e10     # cm s^-1
K_OVER_HC_CM = KB / (H * C_CM_S)  # 0.695034800... cm^-1 K^-1


def parse_state_label(label):
    """Return (v, j) from a label like v1j20."""
    m = re.search(r"v(?P<v>\d+)j(?P<j>\d+)", label)
    if not m:
        raise ValueError(f"Could not parse state label: {label!r}")
    return int(m.group("v")), int(m.group("j"))


def state_label(v, j):
    return f"v{v}j{j}"


def parse_input_filename(path):
    """
    Parse filenames that contain a state label anywhere in the basename.

    Examples accepted:
        v1j0ExpHP.rc
        v1j0_dv0.s
        e1027_v1j0_dv0_avgct.cbt

    Returns a dictionary with prefix, v, j, tail, basename, path.
    """
    base = os.path.basename(path)
    m = re.search(r"(?P<prefix>.*)v(?P<v>\d+)j(?P<j>\d+)(?P<tail>.*)$", base)
    if not m:
        raise ValueError(
            f"Could not find a state label like v1j0 in filename: {base}"
        )
    return {
        "prefix": m.group("prefix"),
        "v": int(m.group("v")),
        "j": int(m.group("j")),
        "tail": m.group("tail"),
        "basename": base,
        "path": os.path.abspath(path),
    }


def expand_input_files(patterns):
    """Expand shell-style glob patterns, preserving explicitly given files."""
    files = []
    for item in patterns:
        matches = sorted(glob.glob(item))
        if matches:
            files.extend(matches)
        elif os.path.exists(item):
            files.append(item)
        else:
            raise FileNotFoundError(f"No file matched input pattern: {item}")

    # De-duplicate while preserving order.
    seen = set()
    unique = []
    for f in files:
        af = os.path.abspath(f)
        if af not in seen:
            seen.add(af)
            unique.append(af)
    return unique


def read_energies(path):
    """
    Read energies from an external text file.

    Accepted forms:
        v1j0   379.98509
        1      0      379.98509

    Energies are returned as {(v, j): energy_value}.
    """
    energies = {}
    with open(path) as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()

            try:
                if re.search(r"v\d+j\d+", parts[0]):
                    v, j = parse_state_label(parts[0])
                    e = float(parts[1])
                else:
                    v = int(parts[0])
                    j = int(parts[1])
                    e = float(parts[2])
            except Exception as exc:
                raise ValueError(
                    f"Could not parse energy line {lineno} in {path}: {raw.rstrip()}"
                ) from exc

            energies[(v, j)] = e
    return energies


def read_cs_file(path):
    """
    Read a cross-section file.

    Returns {j_f: (sigma, err)}. Header/text lines are skipped.
    """
    data = {}
    with open(path) as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            try:
                jf = int(float(parts[0]))
                sigma = float(parts[1])
                err = float(parts[2]) if len(parts) >= 3 else 0.0
            except (ValueError, IndexError):
                continue
            data[jf] = (sigma, err)
    return data


def boltzmann_factor(energy, temperature, energy_unit):
    """Return exp(-E/kBT)."""
    if energy_unit == "cm-1":
        exponent = -energy / (K_OVER_HC_CM * temperature)
    elif energy_unit == "J":
        exponent = -energy / (KB * temperature)
    else:
        raise ValueError(f"Unknown energy unit: {energy_unit}")
    return math.exp(exponent)


def fmt(x):
    if x is None:
        return "N/A"
    if isinstance(x, str):
        return x
    if isinstance(x, int):
        return str(x)
    return f"{x:.10g}"


def safe_ratio(num, den):
    if num is None or den is None or den == 0:
        return None
    return num / den


def detailed_balance_for_pair(v, ji, jf, files_by_state, cs_by_state,
                              energies, temperature, energy_unit):
    """Calculate one i->f / f->i detailed-balance ratio."""
    state_i = (v, ji)
    state_f = (v, jf)

    if state_i not in energies:
        raise KeyError(f"Missing energy for {state_label(v, ji)}")
    if state_f not in energies:
        raise KeyError(f"Missing energy for {state_label(v, jf)}")

    sigma_if, err_if = cs_by_state[state_i][jf]
    sigma_fi, err_fi = cs_by_state[state_f][ji]

    Ei = energies[state_i]
    Ef = energies[state_f]

    gi = 2 * ji + 1
    gf = 2 * jf + 1

    bf_i = boltzmann_factor(Ei, temperature, energy_unit)
    bf_f = boltzmann_factor(Ef, temperature, energy_unit)

    forward_weight = gi * bf_i * sigma_if
    reverse_weight = gf * bf_f * sigma_fi

    ratio = safe_ratio(forward_weight, reverse_weight)

    if ratio is not None and sigma_if != 0 and sigma_fi != 0:
        err_ratio = abs(ratio) * math.sqrt((err_if / sigma_if) ** 2 +
                                           (err_fi / sigma_fi) ** 2)
    else:
        err_ratio = None

    return {
        "v": v,
        "j_i": ji,
        "state_i": state_label(v, ji),
        "file_i": os.path.basename(files_by_state[state_i]),
        "E_i": Ei,
        "g_i": gi,
        "boltz_i": bf_i,
        "sigma_if": sigma_if,
        "err_if": err_if,
        "j_f": jf,
        "state_f": state_label(v, jf),
        "file_f": os.path.basename(files_by_state[state_f]),
        "E_f": Ef,
        "g_f": gf,
        "boltz_f": bf_f,
        "sigma_fi": sigma_fi,
        "err_fi": err_fi,
        "forward_weight": forward_weight,
        "reverse_weight": reverse_weight,
        "DB_ratio": ratio,
        "err_DB_ratio": err_ratio,
        "delta_j": jf - ji,
    }


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k)) for k in fieldnames})


def write_ratio_matrix(path, rows, j_values):
    ratio_by_pair = {(r["j_i"], r["j_f"]): r["DB_ratio"] for r in rows}
    err_by_pair = {(r["j_i"], r["j_f"]): r["err_DB_ratio"] for r in rows}

    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["DB_ratio j_i \\ j_f"] + j_values)
        for ji in j_values:
            w.writerow([ji] + [fmt(ratio_by_pair.get((ji, jf))) for jf in j_values])
        w.writerow([])
        w.writerow(["err_DB_ratio j_i \\ j_f"] + j_values)
        for ji in j_values:
            w.writerow([ji] + [fmt(err_by_pair.get((ji, jf))) for jf in j_values])


def write_summary(path, rows):
    by_delta = defaultdict(list)
    for r in rows:
        if r["DB_ratio"] is not None:
            by_delta[abs(r["delta_j"])].append(r)

    out_rows = []
    for dj in sorted(by_delta):
        vals = [r["DB_ratio"] for r in by_delta[dj]]
        errs = [r["err_DB_ratio"] for r in by_delta[dj] if r["err_DB_ratio"] is not None]
        n = len(vals)
        out_rows.append({
            "delta_j": dj,
            "n_pairs": n,
            "average_ratio": mean(vals),
            "rms_propagated_error": math.sqrt(sum(e * e for e in errs)) / n if errs else None,
            "sample_stdev": stdev(vals) if n > 1 else 0.0,
            "min_ratio": min(vals),
            "max_ratio": max(vals),
        })

    write_tsv(path, out_rows,
              ["delta_j", "n_pairs", "average_ratio", "rms_propagated_error",
               "sample_stdev", "min_ratio", "max_ratio"])


def write_by_initial(output_dir, rows):
    by_initial = defaultdict(list)
    for r in rows:
        by_initial[(r["v"], r["j_i"])].append(r)

    for (v, ji), group in sorted(by_initial.items()):
        path = os.path.join(output_dir, f"DBji{ji}.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["j_f", "DB_ratio", "err_DB_ratio"])
            for r in sorted(group, key=lambda x: x["j_f"]):
                w.writerow([r["j_f"], fmt(r["DB_ratio"]), fmt(r["err_DB_ratio"])])


def main():
    p = argparse.ArgumentParser(
        description="Calculate Excel-style detailed balance ratios for a set of dv=0 rotational files."
    )
    p.add_argument("--energies", required=True,
                   help="External energy file: either 'v1j0 E' or 'v j E'.")
    p.add_argument("--inputs", nargs="+", required=True,
                   help="Input cross-section files or glob patterns, e.g. 'v1j*ExpHP.rc'.")
    p.add_argument("--temperature", "-T", type=float, default=636.06,
                   help="Effective temperature in K. Default: 636.06, matching the spreadsheet.")
    p.add_argument("--energy-unit", choices=["cm-1", "J"], default="cm-1",
                   help="Energy unit in the external energy file. Default: cm-1.")
    p.add_argument("--output-dir", default=".",
                   help="Directory for output files. Default: current directory.")
    p.add_argument("--pair-mode", choices=["upward", "all"], default="upward",
                   help="'upward' uses only j_i < j_f like the spreadsheet; 'all' writes both directions.")
    p.add_argument("--no-by-initial", action="store_true",
                   help="Do not write DBji*.csv files.")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    energies = read_energies(args.energies)
    input_files = expand_input_files(args.inputs)

    files_by_state = {}
    cs_by_state = {}
    for path in input_files:
        info = parse_input_filename(path)
        state = (info["v"], info["j"])
        if state in files_by_state:
            raise ValueError(
                f"Duplicate file for {state_label(*state)}: {files_by_state[state]} and {path}"
            )
        files_by_state[state] = path
        cs_by_state[state] = read_cs_file(path)

    v_values = sorted({v for v, _ in files_by_state})
    if len(v_values) != 1:
        raise ValueError(f"Expected one vibrational state in the input set, found v={v_values}")
    v = v_values[0]

    j_values = sorted(j for vv, j in files_by_state if vv == v)

    rows = []
    warnings = []
    for ji in j_values:
        for jf in j_values:
            if ji == jf:
                continue
            if args.pair_mode == "upward" and ji > jf:
                continue
            if jf not in cs_by_state[(v, ji)]:
                warnings.append(f"Missing sigma({state_label(v, ji)} -> j{jf}) in {os.path.basename(files_by_state[(v, ji)])}")
                continue
            if ji not in cs_by_state[(v, jf)]:
                warnings.append(f"Missing sigma({state_label(v, jf)} -> j{ji}) in {os.path.basename(files_by_state[(v, jf)])}")
                continue
            rows.append(detailed_balance_for_pair(
                v, ji, jf, files_by_state, cs_by_state,
                energies, args.temperature, args.energy_unit
            ))

    pair_fields = [
        "v", "j_i", "state_i", "file_i", "E_i", "g_i", "boltz_i",
        "sigma_if", "err_if",
        "j_f", "state_f", "file_f", "E_f", "g_f", "boltz_f",
        "sigma_fi", "err_fi",
        "forward_weight", "reverse_weight", "DB_ratio", "err_DB_ratio", "delta_j",
    ]

    pairs_path = os.path.join(args.output_dir, "detailed_balance_pairs.tsv")
    matrix_path = os.path.join(args.output_dir, "detailed_balance_matrix.tsv")
    summary_path = os.path.join(args.output_dir, "detailed_balance_summary_by_delta_j.tsv")
    warnings_path = os.path.join(args.output_dir, "detailed_balance_warnings.txt")

    write_tsv(pairs_path, rows, pair_fields)
    write_ratio_matrix(matrix_path, rows, j_values)
    write_summary(summary_path, rows)
    if not args.no_by_initial:
        write_by_initial(args.output_dir, rows)

    with open(warnings_path, "w") as f:
        if warnings:
            f.write("\n".join(warnings) + "\n")
        else:
            f.write("No missing reverse/forward pairs found.\n")

    print(f"Temperature: {args.temperature:g} K")
    print(f"Energy unit: {args.energy_unit}")
    print(f"Input files: {len(input_files)}")
    print(f"Detailed-balance pairs written: {len(rows)}")
    print(f"Wrote: {pairs_path}")
    print(f"Wrote: {matrix_path}")
    print(f"Wrote: {summary_path}")
    if not args.no_by_initial:
        print(f"Wrote: {os.path.join(args.output_dir, 'DBji*.csv')}")
    print(f"Wrote: {warnings_path}")


if __name__ == "__main__":
    main()
