#!/usr/bin/env python3
"""
MRrot_dv0.py

Microreversibility / detailed-balance calculator for rotationally inelastic dv = 0 cases.

Works with filename formats such as:

    v1j0_dv0.s
    v1j2_dv0.s
    e1027_v1j0_dv0.s
    e1027_v1j2_dv0.s
    e1027_v1j4_dv0_avgt.cbt

Example:
    python3 MRrot_dv0.py --energies energies.txt --forward e1027_v1j2_dv0.s

For forward file:

    e1027_v1j2_dv0.s

If the forward rows contain:

    jf = 0, 2, 4, 6, ...

then the code reads reverse cross sections from:

    e1027_v1j0_dv0.s
    e1027_v1j2_dv0.s
    e1027_v1j4_dv0.s
    e1027_v1j6_dv0.s

From each reverse file, it extracts the row where the first column equals
the original forward initial j value.

Example:
    forward file = e1027_v1j2_dv0.s
    original ji  = 2

Then for every reverse file, the code extracts the row:

    jf = 2

Detailed-balance ratio:

    ratio = [(2*ji + 1) * E_i * sigma_if] /
            [(2*jf + 1) * E_f * sigma_fi]

Error propagation:

    err_ratio = |ratio| * sqrt(
        (err_forward / forward)^2 + (err_reverse / reverse)^2
    )
"""

import os
import re
import argparse
import math


def read_energies(path):
    """
    Reads a two-column energy file.

    Expected format:

        v1j0    123.456
        v1j2    124.789
        v1j4    127.000

    Blank lines, comments, and non-data lines are skipped.
    """
    energies = {}

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 2:
                continue

            state = parts[0]

            try:
                energy = float(parts[1])
            except ValueError:
                continue

            energies[state] = energy

    return energies


def parse_forward_filename(filename):
    """
    Parses forward filenames such as:

        v1j0_dv0.s
        v1j2_dv0.s
        e1027_v1j0_dv0.s
        e1027_v1j2_dv0.s
        e1027_v1j4_dv0_avgt.cbt

    Important part:

        v<vi>j<ji>_dv0

    Extra prefix before v... is preserved.
    Extra suffix after dv0 is preserved.

    Example:

        e1027_v1j2_dv0.s

    gives:

        prefix = e1027_
        vi     = 1
        ji     = 2
        dv     = 0
        suffix = .s
    """
    base = os.path.basename(filename)

    m = re.search(
        r"(?P<prefix>.*)v(?P<vi>\d+)j(?P<ji>\d+)_dv(?P<dv>-?\d+)(?P<suffix>.*)$",
        base,
    )

    if not m:
        raise ValueError(
            f"Could not parse filename '{base}'. "
            "Expected something like v1j0_dv0.s or e1027_v1j0_dv0.s"
        )

    prefix = m.group("prefix")
    vi = int(m.group("vi"))
    ji = int(m.group("ji"))
    dv = int(m.group("dv"))
    suffix = m.group("suffix")

    if dv != 0:
        raise ValueError(
            f"This script is only for dv=0, but the input file has dv={dv}"
        )

    return prefix, vi, ji, dv, suffix


def read_cs_file(path, target_j=None):
    """
    Reads a cross-section file.

    Accepted data rows:

        j_f    sigma
        j_f    sigma    error

    Header lines, comments, blank lines, and non-data lines are skipped.

    If target_j is None:
        returns all rows as:

            [(jf, sigma, err), ...]

    If target_j is given:
        returns only:

            sigma, err

        for the matching jf row.
    """
    rows = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            try:
                jf = int(float(parts[0]))
                sigma = float(parts[1])
                err = float(parts[2]) if len(parts) >= 3 else 0.0
            except (ValueError, IndexError):
                continue

            if target_j is not None and jf == target_j:
                return sigma, err

            rows.append((jf, sigma, err))

    if target_j is not None:
        raise ValueError(f"jf={target_j} was not found in {path}")

    return rows


def find_reverse_file(base_dir, prefix, vi, jf, suffix):
    """
    Builds the reverse filename using the same prefix and suffix as the forward file.

    Example:

        forward file:
            e1027_v1j2_dv0.s

        forward row:
            jf = 4

        reverse file searched:
            e1027_v1j4_dv0.s

    If suffix is '_avgt.cbt', then it searches:

            e1027_v1j4_dv0_avgt.cbt
    """
    reverse_name = f"{prefix}v{vi}j{jf}_dv0{suffix}"
    reverse_path = os.path.join(base_dir, reverse_name)

    if os.path.exists(reverse_path):
        return reverse_path

    raise FileNotFoundError(f"Missing reverse file: {reverse_name}")


def make_reverse_output_name(forward_file):
    """
    Creates a reverse-output filename.

    Example:

        e1027_v1j2_dv0.s

    becomes:

        e1027_v1j2_dv0_rev.s
    """
    base_dir = os.path.dirname(os.path.abspath(forward_file))
    base = os.path.basename(forward_file)
    root, ext = os.path.splitext(base)

    return os.path.join(base_dir, f"{root}_rev{ext}")


def make_ratio_output_name(forward_file, ji):
    """
    Creates the ratio-output filename.

    Example:

        forward ji = 2

    gives:

        MRji2.csv
    """
    base_dir = os.path.dirname(os.path.abspath(forward_file))

    return os.path.join(base_dir, f"MRji{ji}.csv")


def fmt(x):
    """
    Formats numbers for output.

    None becomes N/A.
    """
    if x is None:
        return "N/A"

    return f"{x:.6f}"


def safe_ratio_error(ratio, forward_val, err_forward, reverse_val, err_reverse):
    """
    Propagates uncertainty for:

        ratio = forward / reverse

    where:

        forward = (2ji + 1) Ei sigma_if
        reverse = (2jf + 1) Ef sigma_fi
    """
    if ratio is None:
        return None

    if forward_val == 0 or reverse_val == 0:
        return None

    return abs(ratio) * math.sqrt(
        (err_forward / forward_val) ** 2
        + (err_reverse / reverse_val) ** 2
    )


def main(energies_file, forward_file):
    energies = read_energies(energies_file)

    prefix, vi, ji, dv, suffix = parse_forward_filename(forward_file)

    state_i = f"v{vi}j{ji}"

    if state_i not in energies:
        raise KeyError(
            f"Energy for initial state {state_i} was not found in {energies_file}"
        )

    energy_i = energies[state_i]
    stat_i = 2 * ji + 1

    forward_rows = read_cs_file(forward_file)

    base_dir = os.path.dirname(os.path.abspath(forward_file)) or os.getcwd()

    main_output = os.path.join(base_dir, "results.tsv")
    reverse_output = make_reverse_output_name(forward_file)
    ratio_output = make_ratio_output_name(forward_file, ji)

    with open(main_output, "w") as out, \
         open(reverse_output, "w") as rev_out, \
         open(ratio_output, "w") as ratio_out:

        out.write(
            "j_i\tstate_i\tE_i\tsigma_if\terr_if\t"
            "j_f\tstate_f\tE_f\treverse_file\treverse_target_j\t"
            "sigma_fi\terr_fi\tforward\treverse\tratio\terr_ratio\tstatus\n"
        )

        rev_out.write(
            "j_f\tsigma_fi\terr_fi\treverse_file\treverse_target_j\tstatus\n"
        )

        ratio_out.write(
            "j_f\tratio\terr_ratio\tstatus\n"
        )

        for jf, sigma_if, err_if in forward_rows:
            state_f = f"v{vi}j{jf}"
            energy_f = energies.get(state_f)

            forward_val = stat_i * energy_i * sigma_if
            err_forward = stat_i * energy_i * err_if

            reverse_file_used = None
            reverse_target_j = ji

            sigma_fi = None
            err_fi = None

            reverse_val = None
            err_reverse = None

            ratio = None
            err_ratio = None

            status = "OK"

            if energy_f is None:
                status = f"Missing energy for {state_f}"

            else:
                stat_f = 2 * jf + 1

                try:
                    reverse_path = find_reverse_file(
                        base_dir,
                        prefix,
                        vi,
                        jf,
                        suffix,
                    )

                    reverse_file_used = os.path.basename(reverse_path)

                    sigma_fi, err_fi = read_cs_file(
                        reverse_path,
                        target_j=reverse_target_j,
                    )

                    reverse_val = stat_f * energy_f * sigma_fi
                    err_reverse = stat_f * energy_f * err_fi

                    if reverse_val == 0:
                        status = "Reverse value is zero"

                    elif forward_val == 0:
                        status = "Forward value is zero"

                    else:
                        ratio = forward_val / reverse_val

                        err_ratio = safe_ratio_error(
                            ratio,
                            forward_val,
                            err_forward,
                            reverse_val,
                            err_reverse,
                        )

                except FileNotFoundError as e:
                    status = str(e)

                except ValueError as e:
                    status = str(e)

            out.write(
                "\t".join(
                    [
                        str(ji),
                        state_i,
                        fmt(energy_i),
                        fmt(sigma_if),
                        fmt(err_if),
                        str(jf),
                        state_f,
                        fmt(energy_f),
                        reverse_file_used if reverse_file_used is not None else "N/A",
                        str(reverse_target_j),
                        fmt(sigma_fi),
                        fmt(err_fi),
                        fmt(forward_val),
                        fmt(reverse_val),
                        fmt(ratio),
                        fmt(err_ratio),
                        status,
                    ]
                )
                + "\n"
            )

            rev_out.write(
                "\t".join(
                    [
                        str(jf),
                        fmt(sigma_fi),
                        fmt(err_fi),
                        reverse_file_used if reverse_file_used is not None else "N/A",
                        str(reverse_target_j),
                        status,
                    ]
                )
                + "\n"
            )

            ratio_out.write(
                "\t".join(
                    [
                        str(jf),
                        fmt(ratio),
                        fmt(err_ratio),
                        status,
                    ]
                )
                + "\n"
            )

    print(f"Forward file: {os.path.basename(forward_file)}")
    print(f"Filename prefix: '{prefix}'")
    print(f"Forward initial state: {state_i}")
    print(f"Forward initial j_i: {ji}")
    print(f"Reverse target row in each reverse file: jf = {ji}")
    print(f"Wrote: {main_output}")
    print(f"Wrote: {reverse_output}")
    print(f"Wrote: {ratio_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="dv=0 rotational microreversibility / detailed-balance calculator."
    )

    parser.add_argument(
        "--energies",
        required=True,
        help="Two-column energy file, e.g. energies.txt",
    )

    parser.add_argument(
        "--forward",
        required=True,
        help="Forward cross-section file, e.g. e1027_v1j2_dv0.s",
    )

    args = parser.parse_args()

    main(args.energies, args.forward)
