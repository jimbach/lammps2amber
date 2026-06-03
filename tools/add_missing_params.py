#!/usr/bin/python3

import re
import os
import argparse
import sys

TORSION_PATTERN = re.compile(r"No torsion terms for\s+(.+)$")
VALID_TORSION = re.compile(r"^[A-Za-z0-9+_-]+(-[A-Za-z0-9+_-]+){3}$")

def fail(msg):
    print(f"[ERROR] {msg}")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Insert missing torsions into frcmod DIHE section"
    )
    parser.add_argument("input_file", help="File containing torsion warnings")
    parser.add_argument("frcmod_file", help="Input frcmod file")
    parser.add_argument("-o", "--output", help="Output frcmod file (default: overwrite input frcmod)")
    args = parser.parse_args()

    input_file = args.input_file
    frcmod_file = args.frcmod_file
    output_file = args.output or frcmod_file

    # ---- file checks ----
    if not os.path.isfile(input_file):
        fail(f"Input file not found: {input_file}")

    if not os.path.isfile(frcmod_file):
        fail(f"Frcmod file not found: {frcmod_file}")

    if args.output and os.path.exists(output_file):
        print(f"[WARNING] Output file exists and will be overwritten: {output_file}")

    # ---- read torsions ----
    with open(input_file, "r") as f:
        lines = f.readlines()

    torsions = []
    for line in lines:
        m = TORSION_PATTERN.search(line)
        if m:
            torsion = "-".join(m.group(1).split())

            if not VALID_TORSION.match(torsion):
                fail(f"Invalid torsion format detected: {torsion}")

            torsions.append(torsion)

    if not torsions:
        fail("No torsions found in input file")

    # ---- read frcmod ----
    with open(frcmod_file, "r") as f:
        frcmod_lines = f.readlines()

    dihe_index = None
    for i, line in enumerate(frcmod_lines):
        if line.strip().upper() == "DIHE":
            dihe_index = i
            break

    if dihe_index is None:
        fail("DIHE section not found in frcmod file")

    # ---- format output ----
    formatted = [
        f"{t:<18}0.0000   0.0000   0.0000\n"
        for t in torsions
    ]

    new_lines = (
        frcmod_lines[:dihe_index + 1]
        + formatted
        + frcmod_lines[dihe_index + 1:]
    )

    # ---- write output ----
    with open(output_file, "w") as f:
        f.writelines(new_lines)

    print(f"[OK] Inserted {len(torsions)} torsions into {output_file}")

if __name__ == "__main__":
    main()
