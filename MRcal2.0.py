#!/usr/bin/env python3
import os
import re
import argparse
import math


def read_energies(path):
    energies = {}

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            state, val = line.split()[:2]
            energies[state] = float(val)

    return energies


def parse_forward_filename(filename):
    base = os.path.basename(filename)

    m = re.match(
        r"^v(?P<vi>\d+)j(?P<ji>\d+)_dv(?P<dv>-?\d+)(?P<suffix>.*)$",
        base
    )

    if not m:
        raise ValueError(
            f"Forward filename must start like v<vi>j<ji>_dv<dv>, got '{base}'"
        )

    return (
        int(m.group("vi")),
        int(m.group("ji")),
        int(m.group("dv")),
        m.group("suffix")
    )


def read_forward_cs(path):
    rows = []

    with open(path) as f:
        for _ in range(5):
            next(f)

        for line in f:
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            jf = int(parts[0])
            sigma = float(parts[1])
            err = float(parts[2])

            rows.append((jf, sigma, err))

    return rows


def read_reverse_cs(state_f, dv_rev, suffix, base_dir, target_j):
    fname = f"{state_f}_dv{dv_rev}{suffix}"
    full = os.path.join(base_dir, fname)

    if not os.path.exists(full):
        raise FileNotFoundError(f"Missing reverse file: {fname}")

    with open(full) as f:
        for _ in range(5):
            next(f)

        for line in f:
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            jf_line = int(parts[0])

            if jf_line == target_j:
                sigma = float(parts[1])
                err = float(parts[2])
                return sigma, err

    raise ValueError(f"jf={target_j} not found in {fname}")


def make_reverse_output_name(forward_file):
    base_dir = os.path.dirname(os.path.abspath(forward_file))
    base = os.path.basename(forward_file)

    root, ext = os.path.splitext(base)
    return os.path.join(base_dir, f"{root}_rev{ext}")


def make_ratio_output_name(forward_file, ji):
    base_dir = os.path.dirname(os.path.abspath(forward_file))
    return os.path.join(base_dir, f"MRji{ji}.csv")


def fmt(x):
    return f"{x:.6f}" if x is not None else "N/A"


def main(energies_file, forward_file):
    energies = read_energies(energies_file)

    vi, ji, dv, suffix = parse_forward_filename(forward_file)

    state_i = f"v{vi}j{ji}"

    if state_i not in energies:
        raise KeyError(f"DeltaE for initial state {state_i} not found")

    delta_E_i = energies[state_i]
    stat_i = 2 * ji + 1

    forward_rows = read_forward_cs(forward_file)

    base_dir = os.path.dirname(os.path.abspath(forward_file)) or os.getcwd()
    dv_rev = -dv

    # AUTO FILE NAMES
    main_output = os.path.join(base_dir, "results.tsv")
    reverse_output = make_reverse_output_name(forward_file)
    ratio_output = make_ratio_output_name(forward_file, ji)

    with open(main_output, "w") as out, \
         open(reverse_output, "w") as rev_out, \
         open(ratio_output, "w") as ratio_out:

        out.write("j_i\tDeltaE_i\tsigma_if\terr_if\tj_f\tDeltaE_f\tsigma_fi\terr_fi\tforward\treverse\tratio\terr_ratio\n")
        rev_out.write("j_f\tsigma_fi\terr_fi\n")
        ratio_out.write("j_f\tratio\terr_ratio\n")

        for jf, sigma_if, err_if in forward_rows:
            forward_val = stat_i * delta_E_i * sigma_if
            err_forward = stat_i * delta_E_i * err_if

            state_f = f"v{vi + dv}j{jf}"
            delta_E_f = energies.get(state_f)

            sigma_fi = err_fi = reverse_val = err_reverse = None
            ratio = err_ratio = None

            if delta_E_f is not None:
                stat_f = 2 * jf + 1

                try:
                    sigma_fi, err_fi = read_reverse_cs(
                        state_f,
                        dv_rev,
                        suffix,
                        base_dir,
                        ji
                    )

                    reverse_val = stat_f * delta_E_f * sigma_fi
                    err_reverse = stat_f * delta_E_f * err_fi

                except:
                    pass

            if (
                reverse_val is not None and
                reverse_val != 0 and
                forward_val != 0
            ):
                ratio = forward_val / reverse_val

                err_ratio = abs(ratio) * math.sqrt(
                    (err_forward / forward_val) ** 2 +
                    (err_reverse / reverse_val) ** 2
                )

            # MAIN FILE
            out.write("\t".join([
                str(ji),
                fmt(delta_E_i),
                fmt(sigma_if),
                fmt(err_if),
                str(jf),
                fmt(delta_E_f),
                fmt(sigma_fi),
                fmt(err_fi),
                fmt(forward_val),
                fmt(reverse_val),
                fmt(ratio),
                fmt(err_ratio)
            ]) + "\n")

            # REVERSE FILE
            rev_out.write("\t".join([
                str(jf),
                fmt(sigma_fi),
                fmt(err_fi)
            ]) + "\n")

            # RATIO FILE
            ratio_out.write("\t".join([
                str(jf),
                fmt(ratio),
                fmt(err_ratio)
            ]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Microreversibility calculator (auto output files)"
    )

    parser.add_argument("--energies", required=True)
    parser.add_argument("--forward", required=True)

    args = parser.parse_args()

    main(args.energies, args.forward)
