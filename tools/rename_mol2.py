#!/usr/bin/env python3

import sys

# ============================================================
# Read ONLY Masses section (safe + strict)
# ============================================================

def read_mapping(filename):

    typeid_to_name = {}
    in_masses = False

    with open(filename) as f:

        for line in f:

            stripped = line.strip()

            if stripped == "Masses":
                in_masses = True
                continue

            if in_masses and stripped.startswith(("Bond", "Angle", "Dihedral",
                                                  "Improper", "Pair", "Atoms")):
                break

            if not in_masses:
                continue

            if "#" not in stripped:
                continue

            left, right = stripped.split("#", 1)

            fields = left.split()
            if len(fields) < 2:
                continue

            typeid = int(fields[0])
            atomname = right.strip().split()[0]

            typeid_to_name[typeid] = atomname

    return typeid_to_name


# ============================================================
# Read atom ID -> atom type name
# ============================================================

def read_atoms(filename, typeid_to_name):

    atomid_to_name = {}
    in_atoms = False

    with open(filename) as f:

        for line in f:

            stripped = line.strip()

            if stripped == "Atoms":
                in_atoms = True
                continue
                
            if stripped == "Bonds":
                in_atoms = False
                continue

            if not in_atoms:
                continue

            if not stripped:
                continue

            fields = stripped.split()

            if len(fields) < 4:
                continue

            atomid = int(fields[0])
            typeid = int(fields[2])

            if typeid in typeid_to_name:
                atomid_to_name[atomid] = typeid_to_name[typeid]

    return atomid_to_name


# ============================================================
# Shorten atom type name
# ============================================================

def shorten(name):
    name = name.replace("_", "").replace(".", "")
    return name[:2]


# ============================================================
# Modify MOL2 (FIXED STRICT FORMAT)
# ============================================================

def modify_mol2(mol2_in, mol2_out, atomid_to_name):

    out = []
    in_atom = False

    with open(mol2_in) as f:

        for line in f:

            s = line.strip()

            if s.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                out.append(line.rstrip("\n"))
                continue

            if s.startswith("@<TRIPOS>BOND"):
                in_atom = False

            if in_atom and s:

                fields = line.split()

                if len(fields) >= 9:

                    atomid = int(fields[0])

                    if atomid in atomid_to_name:

                        newtype = shorten(atomid_to_name[atomid])

                        atom_name = fields[1]

                        x = float(fields[2])
                        y = float(fields[3])
                        z = float(fields[4])

                        subst_id = fields[6]
                        subst_name = fields[7]
                        charge = float(fields[8])

                        # =====================================================
                        # FIXED MOL2 FORMAT (prevents all column merging)
                        # =====================================================
                        newline = (
                            f"{atomid:7d} "        # ID + gap
                            f"{atom_name:4s}     "     # atom name + gap
                            f"{x:10.4f} "            # X
                            f"{y:10.4f} "            # Y
                            f"{z:10.4f} "           # Z + gap
                            f"{newtype:6s} "        # atom type
                            f" {subst_id:>4s} "       # residue id
                            f"{subst_name:6s}    "     # residue name
                            f"{charge:8.6f}"        # charge
                        )

                        out.append(newline)
                        continue

            out.append(line.rstrip("\n"))

    with open(mol2_out, "w") as f:
        f.write("\n".join(out) + "\n")


# ============================================================
# MAIN
# ============================================================

def main():

    if len(sys.argv) != 4:
        print("Usage: python script.py file1.txt input.mol2 output.mol2")
        sys.exit(1)

    file1 = sys.argv[1]
    mol2_in = sys.argv[2]
    mol2_out = sys.argv[3]

    type_map = read_mapping(file1)
    print(type_map)
    atom_map = read_atoms(file1, type_map)
    
    modify_mol2(mol2_in, mol2_out, atom_map)

    print("Done ->", mol2_out)


if __name__ == "__main__":
    main()
