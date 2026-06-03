#!/usr/bin/env python3

import sys


def parse_new_bonds(file2):
    """
    Parse bonds from file 2.

    Expected format:

      Bonds

           1        1        1        4
           2        2        1      557

      Angles

    Columns:
        bond_id  bond_type  atom1  atom2

    bond_type is ignored.
    """

    bonds = []

    with open(file2, "r") as f:
        lines = f.readlines()

    in_bonds = False

    for line in lines:
        stripped = line.strip()

        # Start reading after "Bonds"
        if stripped == "Bonds":
            in_bonds = True
            continue

        if not in_bonds:
            continue

        # Stop reading at "Angles"
        if stripped == "Angles":
            break

        # Skip empty lines
        if not stripped:
            continue

        parts = stripped.split()

        # Skip malformed lines
        if len(parts) < 4:
            continue

        try:
            bond_id = int(parts[0])
            atom1 = int(parts[2])
            atom2 = int(parts[3])
        except ValueError:
            continue

        bonds.append((bond_id, atom1, atom2))

    return bonds


def update_molecule_counts(lines, new_bond_count):
    """
    Update the bond count in the MOLECULE section.

    Example:

      276   300     1     0     0

    becomes:

      276   725     1     0     0
    """

    output = []

    in_molecule_section = False
    updated = False

    for i, line in enumerate(lines):

        stripped = line.strip()

        output.append(line)

        if stripped == "@<TRIPOS>MOLECULE":
            in_molecule_section = True
            continue

        if in_molecule_section and not updated:

            # The counts line is typically 2 lines below @<TRIPOS>MOLECULE
            # and contains 5 integers.
            parts = line.split()

            if len(parts) == 5 and all(p.lstrip("-").isdigit() for p in parts):

                num_atoms = int(parts[0])
                num_subst = int(parts[2])
                num_feat = int(parts[3])
                num_sets = int(parts[4])

                new_line = (
                    f"{num_atoms:5d}"
                    f"{new_bond_count:6d}"
                    f"{num_subst:6d}"
                    f"{num_feat:6d}"
                    f"{num_sets:6d}\n"
                )

                output[-1] = new_line

                updated = True
                in_molecule_section = False

    if not updated:
        raise RuntimeError("Could not find MOLECULE counts line.")

    return output


def replace_bond_section(lines, bonds):
    """
    Replace the @<TRIPOS>BOND section while preserving
    all other formatting/content.
    """

    output = []

    in_bond_section = False
    inserted = False

    for line in lines:

        stripped = line.strip()

        # Start of bond section
        if stripped == "@<TRIPOS>BOND":
            in_bond_section = True
            output.append(line)

            # Insert new bonds
            for bond_id, atom1, atom2 in bonds:
                output.append(
                    f"{bond_id:6d}{atom1:6d}{atom2:6d} 1   \n"
                )

            inserted = True
            continue

        # Skip original bond lines
        if in_bond_section:

            # Next mol2 section starts
            if stripped.startswith("@<TRIPOS>") and stripped != "@<TRIPOS>BOND":
                in_bond_section = False
                output.append(line)

            continue

        output.append(line)

    if not inserted:
        raise RuntimeError("No @<TRIPOS>BOND section found.")

    return output


def main():

    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input.mol2 bonds.txt output.mol2")
        sys.exit(1)

    mol2_file = sys.argv[1]
    bonds_file = sys.argv[2]
    output_file = sys.argv[3]

    bonds = parse_new_bonds(bonds_file)

    if not bonds:
        print("Error: No bonds found in second file.")
        sys.exit(1)

    with open(mol2_file, "r") as f:
        lines = f.readlines()

    # Update bond count in MOLECULE section
    lines = update_molecule_counts(lines, len(bonds))

    # Replace BOND section
    lines = replace_bond_section(lines, bonds)

    with open(output_file, "w") as f:
        f.writelines(lines)

    print(f"Wrote updated mol2 file: {output_file}")


if __name__ == "__main__":
    main()
