#!/usr/bin/env python3
import os
import re
import argparse
import math


def read_energies(path):
    """Load 'state  ΔE' lines into a dict (state → ΔE)."""
    energies = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            state, val = line.split()
            energies[state] = float(val)
    return energies


def parse_forward_filename(filename):
    """
    Parse a forward file named like 'v1j0.-1s'.
    Extracts:
      vi: vibrational index (int)
      ji: rotational index (int)
      dv: vibrational change (int)
      method: binning method ('s' or 'n')
    Returns (vi, ji, dv, method).
    """
    base = os.path.basename(filename)
    m = re.match(r"^v(?P<vi>\d+)j(?P<ji>\d+)\.(?P<dv>-?\d+)(?P<method>[sn])$", base)
    if not m:
        raise ValueError(
            f"Forward filename must be v<vi>j<ji>.<dv><method>, got '{base}'"
        )
    return int(m.group('vi')), int(m.group('ji')), int(m.group('dv')), m.group('method')


def read_forward_cs(path):
    """
    Reads forward cross-section file, skipping first five header lines,
    then data lines 'jf sigma error'.
    Returns list of (jf, sigma, err)."""
    rows = []
    with open(path) as f:
        for _ in range(5):
            next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            jf, sigma, err = int(parts[0]), float(parts[1]), float(parts[2])
            rows.append((jf, sigma, err))
    return rows


def read_reverse_cs(state_f, dv_rev, method, base_dir, target_j):
    """
    Reads reverse cross-section file '<state_f>.<dv_rev><method>', skipping first five lines.
    Finds the line where jf == target_j, returns (sigma, err).
    """
    fname = f"{state_f}.{dv_rev}{method}"
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
                err   = float(parts[2])
                return sigma, err
    raise ValueError(f"jf={target_j} not found in {fname}")


def main(energies_file, forward_file, output_file):
    # Load ΔE values
    energies = read_energies(energies_file)

    # Parse forward file name
    vi, ji, dv, method = parse_forward_filename(forward_file)
    state_i = f"v{vi}j{ji}"
    if state_i not in energies:
        raise KeyError(f"ΔE for initial state {state_i} not found in {energies_file}")
    ΔE_i = energies[state_i]
    stat_i = 2 * ji + 1

    # Read forward data: jf, sigma_if, err_if
    forward_rows = read_forward_cs(forward_file)

    # Directory for reverse files
    base_dir = os.path.dirname(os.path.abspath(forward_file)) or os.getcwd()
    dv_rev = -dv

    # Write header with extra columns
    with open(output_file, 'w') as out:
        header = [
            'j_i', 'ΔE_i', 'sigma_if', 'err_if',
            'j_f', 'ΔE_f', 'sigma_fi', 'err_fi',
            'forward', 'reverse', 'ratio', 'err_ratio'
        ]
        out.write('\t'.join(header) + '\n')

        # Process each j_f
        for jf, sigma_if, err_if in forward_rows:
            # Forward component and error
            forward_val = stat_i * ΔE_i * sigma_if
            err_forward = stat_i * ΔE_i * err_if

            # Reverse component and error for transition back to state_i
            state_f = f"v{vi+dv}j{jf}"
            ΔE_f = energies.get(state_f)
            if ΔE_f is None:
                sigma_fi = err_fi = reverse_val = err_reverse = None
            else:
                stat_f = 2 * jf + 1
                try:
                    sigma_fi, err_fi = read_reverse_cs(state_f, dv_rev, method, base_dir, ji)
                    reverse_val = stat_f * ΔE_f * sigma_fi
                    err_reverse = stat_f * ΔE_f * err_fi
                except (FileNotFoundError, ValueError):
                    sigma_fi = err_fi = reverse_val = err_reverse = None

            # Ratio and error propagation
            if reverse_val and reverse_val != 0:
                ratio = forward_val / reverse_val
                err_ratio = abs(ratio) * math.sqrt(
                    (err_forward / forward_val) ** 2 +
                    (err_reverse / reverse_val) ** 2
                )
            else:
                ratio = None
                err_ratio = None

            # Format values or N/A
            def fmt(x): return f"{x:.6f}" if x is not None else 'N/A'

            row = [
                str(ji), fmt(ΔE_i), fmt(sigma_if), fmt(err_if),
                str(jf), fmt(ΔE_f), fmt(sigma_fi), fmt(err_fi),
                fmt(forward_val), fmt(reverse_val), fmt(ratio), fmt(err_ratio)
            ]
            out.write('\t'.join(row) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compute forward, reverse, ratio, and ratio error to a TSV file"
    )
    parser.add_argument('--energies', required=True,
                        help="Two-column ΔE file, e.g. energies.txt")
    parser.add_argument('--forward', required=True,
                        help="Forward file named v<vi>j<ji>.<dv><method>, e.g. v1j2.-1s")
    parser.add_argument('--output', '-o', required=True,
                        help="Output file to write tab-separated results")
    args = parser.parse_args()
    main(args.energies, args.forward, args.output)

