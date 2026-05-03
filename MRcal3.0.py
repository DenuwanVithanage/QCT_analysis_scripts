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
    """
    Supports:
      v1j0_dv-1.s
      v1j0_dv-1_.s
      v1j0_dv-1_avgt.cbt
      v1j0_dv-1avgct.cbt
      v0j10_dv1.mol
    """
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


def read_cs_file(path, target_j=None):
    """
    Flexible cross-section reader.

    Accepts:
      j_f   sigma
      j_f   sigma   error

    Skips headers/comments/non-data lines automatically.
    If error is missing, err = 0.0.
    """
    rows = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            try:
                jf = int(parts[0])
                sigma = float(parts[1])
                err = float(parts[2]) if len(parts) >= 3 else 0.0
            except (ValueError, IndexError):
                continue

            if target_j is not None and jf == target_j:
                return sigma, err

            rows.append((jf, sigma, err))

    if target_j is not None:
        raise ValueError(f"jf={target_j} not found in {path}")

    return rows


def read_reverse_cs(state_f, dv_rev, suffix, base_dir, target_j):
    fname = f"{state_f}_dv{dv_rev}{suffix}"
    full = os.path.join(base_dir, fname)

    if not os.path.exists(full):
        raise FileNotFoundError(f"Missing reverse file: {fname}")

    return read_cs_file(full, target_j=target_j)


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
        raise KeyError(f"DeltaE for initial state {state_i} not found in {energies_file}")

    delta_E_i = energies[state_i]
    stat_i = 2 * ji + 1

    forward_rows = read_cs_file(forward_file)

    base_dir = os.path.dirname(os.path.abspath(forward_file)) or os.getcwd()
    dv_rev = -dv

    main_output = os.path.join(base_dir, "results.tsv")
    reverse_output = make_reverse_output_name(forward_file)
    ratio_output = make_ratio_output_name(forward_file, ji)

    with open(main_output, "w") as out, \
         open(reverse_output, "w") as rev_out, \
         open(ratio_output, "w") as ratio_out:

        out.write(
            "j_i\tDeltaE_i\tsigma_if\terr_if\t"
            "j_f\tDeltaE_f\tsigma_fi\terr_fi\t"
            "forward\treverse\tratio\terr_ratio\n"
        )

        rev_out.write("j_f\tsigma_fi\terr_fi\n")
        ratio_out.write("j_f\tratio\terr_ratio\n")

        for jf, sigma_if, err_if in forward_rows:
            forward_val = stat_i * delta_E_i * sigma_if
            err_forward = stat_i * delta_E_i * err_if

            state_f = f"v{vi + dv}j{jf}"
            delta_E_f = energies.get(state_f)

            sigma_fi = None
            err_fi = None
            reverse_val = None
            err_reverse = None
            ratio = None
            err_ratio = None

            if delta_E_f is not None:
                stat_f = 2 * jf + 1

                try:
                    sigma_fi, err_fi = read_reverse_cs(
                        state_f=state_f,
                        dv_rev=dv_rev,
                        suffix=suffix,
                        base_dir=base_dir,
                        target_j=ji
                    )

                    reverse_val = stat_f * delta_E_f * sigma_fi
                    err_reverse = stat_f * delta_E_f * err_fi

                except (FileNotFoundError, ValueError):
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
                fmt(err_ratio),
            ]) + "\n")

            rev_out.write("\t".join([
                str(jf),
                fmt(sigma_fi),
                fmt(err_fi),
            ]) + "\n")

            ratio_out.write("\t".join([
                str(jf),
                fmt(ratio),
                fmt(err_ratio),
            ]) + "\n")

    print(f"Wrote: {main_output}")
    print(f"Wrote: {reverse_output}")
    print(f"Wrote: {ratio_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Microreversibility calculator with flexible file parsing."
    )

    parser.add_argument(
        "--energies",
        required=True,
        help="Two-column energy file, e.g. energies.txt"
    )

    parser.add_argument(
        "--forward",
        required=True,
        help="Forward file, e.g. v1j0_dv-1.s or v0j10_dv1.mol"
    )

    args = parser.parse_args()

    main(args.energies, args.forward)
