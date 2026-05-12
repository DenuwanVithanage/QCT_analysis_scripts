#!/usr/bin/env python3
import glob, re

# how many header lines (comments/info) to skip in each file
HEADER_LINES = 5

# the targets (column‑1 values) you want to scan
TARGETS    = list(range(0, 61, 2))   # 0,2,4,...,44

# the file‐extensions to process
EXTENSIONS = [".s", "_avgct.cbt", "_avgt.cbt", "_opt.cbt", ".sgut"]

def find_cs_and_error(fn, target, header_lines=HEADER_LINES):
    """
    Open fn, skip header_lines lines, then find a row whose first column == target.
    Return (cross_section, error) from columns 2 and 3, or None if no match.
    """
    with open(fn) as f:
        # skip header
        for _ in range(header_lines):
            next(f, None)

        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0] == target:
                return float(parts[1]), float(parts[2])
    return None

def process_one(target, ext):
    pattern = f"e*_v1j0_dv-1{ext}"
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[WARN] no files match '{pattern}'")
        return

    out_rows = []
    for fn in files:
        m = re.match(r"e(\d+)_", fn)
        if not m:
            continue
        energy_id = int(m.group(1))
        vals = find_cs_and_error(fn, str(target))
        if vals:
            cs, err = vals
            out_rows.append((energy_id, cs, err))

    if not out_rows:
        print(f"[WARN] target={target}: no data found in any .{ext} files")
        return

    # sort by energy ID
    out_rows.sort(key=lambda r: r[0])

    out_fname = f"CS_vs_E_jf_{target}{ext}"
    with open(out_fname, "w") as out:
        for eid, cs, err in out_rows:
            out.write(f"{eid}\t{cs}\t{err}\n")

    print(f"[ OK ] wrote {len(out_rows)} lines → {out_fname}")

def main():
    for target in TARGETS:
        for ext in EXTENSIONS:
            process_one(target, ext)

if __name__ == "__main__":
    main()

