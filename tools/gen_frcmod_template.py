#!/usr/bin/env python3

import sys
from itertools import product


def parse_atom_types(mol2_lines):
    types = []

    in_atom = False

    for line in mol2_lines:
        if line.strip() == "@<TRIPOS>ATOM":
            in_atom = True
            continue

        if in_atom and line.startswith("@<TRIPOS>"):
            break

        if in_atom:
            parts = line.split()
            if len(parts) < 6:
                continue
            types.append(parts[5])

    return types


def trim(t):
    return t[:2]


# ---------- BONDS (unordered OK) ----------
def bonds(types):
    seen = set()
    out = []

    for a, b in product(types, repeat=2):
        key = tuple(sorted((a, b)))
        if key not in seen:
            seen.add(key)
            out.append((a, b))

    return out


# ---------- ANGLES (ordered but symmetric removal optional) ----------
def angles(types):
    seen = set()
    out = []

    for a, b, c in product(types, repeat=3):
        key = (a, b, c)  # KEEP ORDER but avoid exact duplicates only
        if key not in seen:
            seen.add(key)
            out.append((a, b, c))

    return out


# ---------- DIHEDRALS (ORDERED - IMPORTANT) ----------
def dihedrals(types):
    seen = set()
    out = []

    for a, b, c, d in product(types, repeat=4):
        key = (a, b, c, d)
        if key not in seen:
            seen.add(key)
            out.append((a, b, c, d))

    return out


# ---------- formatting ----------
def fmt2(x):
    return trim(x)


def fmt_pair(a, b):
    return f"{fmt2(a):<2s}-{fmt2(b):<2s}"


def fmt_triplet(a, b, c):
    return f"{fmt2(a):<2s}-{fmt2(b):<2s}-{fmt2(c):<2s}"


def fmt_quad(a, b, c, d):
    return f"{fmt2(a):<2s}-{fmt2(b):<2s}-{fmt2(c):<2s}-{fmt2(d):<2s}"


def main():

    if len(sys.argv) != 3:
        print("Usage: script mol2 input.mol2 output.frcmod")
        sys.exit(1)

    mol2_file, out_file = sys.argv[1:3]

    with open(mol2_file) as f:
        lines = f.readlines()

    types = sorted(set(parse_atom_types(lines)))

    if not types:
        raise RuntimeError("No atom types found")

    out = []

    out.append("remark auto-generated frcmod\n\n")

    # MASS
    out.append("MASS\n")
    for t in types:
        out.append(f"{fmt2(t):<4s}   12.0000   0.000\n")
    out.append("\n")

    # BOND
    out.append("BOND\n")
    for a, b in bonds(types):
        out.append(f"{fmt_pair(a,b):<6s}   0.0000   0.0000\n")
    out.append("\n")

    # ANGLE
    out.append("ANGLE\n")
    for a, b, c in angles(types):
        out.append(f"{fmt_triplet(a,b,c):<10s}   0.0000   0.0000\n")
    out.append("\n")

    # DIHE
    out.append("DIHE\n")
    for a, b, c, d in dihedrals(types):
        out.append(f"{fmt_quad(a,b,c,d):<14s}   0.0000   0.0000   0.0000\n")
    out.append("\n")

    # IMPROPER
    out.append("IMPROPER\n")
    for a, b, c, d in dihedrals(types):
        out.append(f"{fmt_quad(a,b,c,d):<14s}   0.0000   0.0000   0.0000\n")
    out.append("\n")

    # NONBON
    out.append("NONBON\n")
    for t in types:
        out.append(f"{fmt2(t):<4s}   0.0000   0.0000\n")

    with open(out_file, "w") as f:
        f.writelines(out)

    print(f"Wrote frcmod: {out_file}")


if __name__ == "__main__":
    main()
