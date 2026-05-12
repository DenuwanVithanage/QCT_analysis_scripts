#!/usr/bin/env python3

import os
import re
import argparse
import math


def read_energies(path):
    """
    Reads a two-column energy file:

        v1j0   1027.123
        v1j2   1020.456
        v1j4   1011.789

    Blank lines and lines starting with # are skipped.
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
            val = float(parts[1])

            energies[state] = val

    return energies


def parse_forward_filename(filename):
    """
    Parses forward filenames such as:

        v1j0_dv0.s
        e1027_v1j0_dv0.s
        e1027_v1j0_dv0_avgt.cbt
        e1027_v1j0_dv0_avgct.cbt

    Returns:

        prefix, vi, ji, dv, suffix

    Example:

        e1027_v1j0_dv0.s

    gives:

        prefix = "e1027_"
        vi     = 1
        ji     = 0
        dv     = 0
        suffix = ".s"
    """
    base = os.path.basename(filename)

    m = re.search(
        r"^(?P<prefix>.*)v(?P<vi>\d+)j(?P<ji>\d+)_dv(?P<dv>-?\d+)(?P<suffix>.*)$",
        base
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
            f"This simplified script is only for rotational dv=0 cases, "
            f"but the input file has dv={dv}"
        )

    return prefix, vi, ji, dv, suffix


def read_cs_file(path, target_j=None):
    """
    Flexible cross-section file reader.

    Accepts files with data like:

        j_f   sigma
        j_f   sigma   error

    It automatically skips:
        - blank lines
        - comment lines beginning with #
        - header/text lines

    If target_j is None:
        returns all rows as:
            [(jf, sigma, err), ...]

    If target_j is given:
        returns:
            sigma, err

    Example:
        read_cs_file("v1j2_dv0.s", target_j=0)

    extracts the reverse cross section for:

        v1j2 -> v1j0
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


def find_reverse_file(base_dir, prefix, vi, jf, suffix):
    """
    For forward file:

        e1027_v1j0_dv0.s

    and forward row:

        jf = 2

    this looks for:

        e1027_v1j2_dv0.s

    The vibrational state does not change because dv = 0.
    """
    reverse_name = f"{prefix}v{vi}j{jf}_dv0{suffix}"
    reverse_path = os.path.join(base_dir, reverse_name)

    if os.path.exists(reverse_path):
        return reverse_path

    raise FileNotFoundError(f"Missing reverse file: {reverse_name}")


def make_reverse_output_name(forward_file):
    """
    Example:

        v1j0_dv0.s          -> v1j0_dv0_rev.s
        e1027_v1j0_dv0.s    -> e1027_v1j0_dv0_rev.s
    """
    base_dir = os.path.dirname(os.path.abspath(forward_file))
    base = os.path.basename(forward_file)

    root, ext = os.path.splitext(base)

    return os.path.join(base_dir, f"{root}_rev{ext}")


def make_ratio_output_name(forward_file, ji):
    """
    Example:

        initial ji = 0  -> MRji0.csv
        initial ji = 4  -> MRji4.csv
    """
    base_dir = os.path.dirname(os.path.abspath(forward_file))
    return os.path.join(base_dir, f"MRji{ji}.csv")


def fmt(x):
    """
    Output formatter.
    """
    if x is None:
        return "N/A"

    return f"{x:.6f}"


def main(energies_file, forward_file):
    energies = read_energies(energies_file)

    prefix, vi, ji, dv, suffix = parse_forward_filename(forward_file)

    state_i = f"v{vi}j{ji}"

    if state_i not in energies:
        raise KeyError(
            f"Energy for initial state {state_i} not found in {energies_file}"
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
            "j_i\tEnergy_i\tsigma_if\terr_if\t"
            "j_f\tEnergy_f\tsigma_fi\terr_fi\t"
            "j_used_forward\tj_used_reverse\t"
            "forward\treverse\tratio\terr_ratio\t"
            "reverse_file\n"
        )

        rev_out.write("j_f\tsigma_fi\terr_fi\n")
        ratio_out.write("j_f\tratio\terr_ratio\n")

        for jf, sigma_if, err_if in forward_rows:

            state_f = f"v{vi}j{jf}"
            energy_f = energies.get(state_f)

            sigma_fi = None
            err_fi = None
            reverse_file_used = "N/A"

            forward_val = None
            reverse_val = None
            err_forward = None
            err_reverse = None
            ratio = None
            err_ratio = None

            if energy_f is not None:

                stat_f = 2 * jf + 1

                forward_val = stat_i * energy_i * sigma_if
                err_forward = stat_i * energy_i * err_if

                try:
                    reverse_path = find_reverse_file(
                        base_dir=base_dir,
                        prefix=prefix,
                        vi=vi,
                        jf=jf,
                        suffix=suffix
                    )

                    reverse_file_used = os.path.basename(reverse_path)

                    # Important:
                    # Forward is v_i j_i -> v_i j_f
                    # Reverse is v_i j_f -> v_i j_i
                    #
                    # Therefore from the reverse file v_i j_f_dv0.s,
                    # we extract the row where final j equals the original ji.
                    sigma_fi, err_fi = read_cs_file(
                        reverse_path,
                        target_j=ji
                    )

                    reverse_val = stat_f * energy_f * sigma_fi
                    err_reverse = stat_f * energy_f * err_fi

                except (FileNotFoundError, ValueError):
                    pass

            if (
                forward_val is not None and
                reverse_val is not None and
                forward_val != 0 and
                reverse_val != 0
            ):
                ratio = forward_val / reverse_val

                err_ratio = abs(ratio) * math.sqrt(
                    (err_forward / forward_val) ** 2 +
                    (err_reverse / reverse_val) ** 2
                )

            out.write("\t".join([
                str(ji),
                fmt(energy_i),
                fmt(sigma_if),
                fmt(err_if),
                str(jf),
                fmt(energy_f),
                fmt(sigma_fi),
                fmt(err_fi),
                str(ji),
                str(jf),
                fmt(forward_val),
                fmt(reverse_val),
                fmt(ratio),
                fmt(err_ratio),
                reverse_file_used,
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
        description="Simplified rotational microreversibility calculator for dv=0."
    )

    parser.add_argument(
        "--energies",
        required=True,
        help="Two-column energy file, e.g. energies.txt"
    )

    parser.add_argument(
        "--forward",
        required=True,
        help="Forward file, e.g. v1j0_dv0.s or e1027_v1j0_dv0.s"
    )

    args = parser.parse_args()

    main(args.energies, args.forward)
