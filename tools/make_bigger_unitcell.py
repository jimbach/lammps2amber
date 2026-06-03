#!/usr/bin/env python3

import sys


def parse_mol2_template(mol2_lines):
    atoms = []
    in_atom = False

    for line in mol2_lines:
        if line.strip() == "@<TRIPOS>ATOM":
            in_atom = True
            continue

        if in_atom and line.startswith("@<TRIPOS>"):
            break

        if in_atom:
            parts = line.split()
            if len(parts) < 9:
                continue

            atoms.append({
                "name": parts[1],
                "type": parts[5],
                "resname": parts[7],
                "charge": float(parts[8]),
            })

    return atoms


def parse_file2(file2_lines):
    coords = []
    in_atoms = False

    for line in file2_lines:
        s = line.strip()

        if s == "Atoms":
            in_atoms = True
            continue

        if not in_atoms:
            continue

        if s.startswith("Bonds"):
            break

        parts = line.split()
        if len(parts) < 7:
            continue

        try:
            x, y, z = map(float, parts[-3:])
        except ValueError:
            continue

        coords.append((x, y, z))

    return coords


def rebuild_atoms(template_atoms, coords):

    out = []

    for i, (x, y, z) in enumerate(coords):
        t = template_atoms[i % len(template_atoms)]

        out.append(
            f"{i+1:7d}"
            f"{t['name']:>7s}"
            f"{x:12.4f}"
            f"{y:12.4f}"
            f"{z:12.4f}"
            f"{t['type']:>6s}"
#            f"{(i // len(template_atoms)) + 1:6d}"
            f"{1:6d}"
            f"{t['resname']:>6s}"
            f"{float(t['charge']):12.6f}\n"
        )

    return out


def update_molecule_block(lines, new_atom_count):

    out = []
    in_mol = False
    updated = False

    for line in lines:

        if line.strip() == "@<TRIPOS>MOLECULE":
            in_mol = True
            out.append(line)
            continue

        if in_mol and not updated:
            parts = line.split()

            # detect counts line (first numeric line after header block starts)
            if len(parts) == 5 and parts[0].isdigit():

                atoms = new_atom_count
                bonds = int(parts[1])
                subst = int(parts[2])
                feat = int(parts[3])
                sets = int(parts[4])

                # fixed-width MOL2 formatting
                new_line = (
                    f"{atoms:5d}"
                    f"{bonds:6d}"
                    f"{subst:6d}"
                    f"{feat:6d}"
                    f"{sets:6d}\n"
                )

                out.append(new_line)
                updated = True
                continue

        out.append(line)

    return out


def replace_atoms(mol2_lines, new_atoms):

    out = []
    in_atom = False

    for line in mol2_lines:

        if line.strip() == "@<TRIPOS>ATOM":
            in_atom = True
            out.append(line)

            for a in new_atoms:
                out.append(a)

            continue

        if in_atom and line.startswith("@<TRIPOS>"):
            in_atom = False
            out.append(line)
            continue

        if in_atom:
            continue

        out.append(line)

    return out


def main():

    if len(sys.argv) != 4:
        print("Usage: script mol2 file2 output.mol2")
        sys.exit(1)

    mol2_file, file2_file, out_file = sys.argv[1:4]

    with open(mol2_file) as f:
        mol2_lines = f.readlines()

    with open(file2_file) as f:
        file2_lines = f.readlines()

    template_atoms = parse_mol2_template(mol2_lines)
    coords = parse_file2(file2_lines)

    print("Template atoms:", len(template_atoms))
    print("File2 atoms:", len(coords))

    new_atoms = rebuild_atoms(template_atoms, coords)

    mol2_lines = replace_atoms(mol2_lines, new_atoms)
    mol2_lines = update_molecule_block(mol2_lines, len(coords))

    with open(out_file, "w") as f:
        f.writelines(mol2_lines)

    print("Done.")


if __name__ == "__main__":
    main()
