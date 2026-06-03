#!/usr/bin/env python3

import re
import sys
import numpy as np
from collections import defaultdict

# ============================================================
# Mapping between LAMMPS atom types and AMBER frcmod atom types
# (UNCHANGED)
# ============================================================

ATOM_MAP = {
    "C_R": "CR",
    "C_3": "C3",
    "H_": "H",
    "N_R": "NR",
    "Zn3f2": "Zn",
}

# ============================================================
# Helpers
# ============================================================

def normalize_key(parts):
    fwd = tuple(parts)
    rev = tuple(reversed(parts))
    return min(fwd, rev)


def map_atoms(atom_list):
    return [ATOM_MAP.get(a, a) for a in atom_list]


# ============================================================
# FORMAT FIX (ONLY CHANGE YOU REQUESTED)
# ============================================================

def fmt_atom(a):
    # ensures H becomes "H "
    return f"{a:<2}"


# ============================================================
# Parse LAMMPS coefficient file
# ============================================================

def parse_lammps_coeffs(filename):

    masses = {}
    bonds = {}
    angles = {}
    dihedrals = {}
    impropers = {}
    nonbond = {}

    section = None

    with open(filename) as f:
        for raw in f:

            line = raw.strip()

            if not line:
                continue

            if line.startswith("Masses"):
                section = "Masses"
                continue
            elif line.startswith("Bond Coeffs"):
                section = "Bond"
                continue
            elif line.startswith("Angle Coeffs"):
                section = "Angle"
                continue
            elif line.startswith("Dihedral Coeffs"):
                section = "Dihedral"
                continue
            elif line.startswith("Improper Coeffs"):
                section = "Improper"
                continue
            elif line.startswith("Pair Coeffs"):
                section = "Pair"
                continue

            if "#" not in line:
                continue

            data, comment = line.split("#", 1)
            atom_types = comment.strip().split()
            fields = data.split()

            # ---------------- MASS ----------------
            if section == "Masses":
                mass = float(fields[1])
                amber = map_atoms([atom_types[0]])[0]
                masses[amber] = mass

            # ---------------- BONDS ----------------
            elif section == "Bond":
                k = float(fields[1])
                req = float(fields[2])

                atoms = map_atoms(atom_types)
                key = normalize_key(atoms)

                bonds[key] = (k, req)

            # ---------------- ANGLES ----------------
            elif section == "Angle":
                style = fields[1]
                atoms = map_atoms(atom_types)
                key = normalize_key(atoms)
                #cosine/periodic 57.289016 -1 3 # H_ C_R C_
                #
                #which means 2/3^2 * 57 (1 - - 1(-1)^3 cos(3 phi))
                #
                #to amber cos^2 like potential (9 from dElammps/dPhi = 3 * 3 *cos(3phi)
                # dElammps/dPhi | 120 deg = 2* 57 / 3^2 * 9 * cos (3phi) = 2* 57 / 3^2 * 9 * cos(3 * 120/3)
                # dElammps/dPhi | 360/n deg = 2* forceconstlammps / n^2 * n^2 
                # dElammps/dPhi | 360/n deg = 2* forceconstlammps 
                # equilibrium at 360 / n 
                
                
                if style == "cosine/periodic":
                    angles[key] = (float(fields[2])*2.0, 360.0/float(fields[4]))

                elif style == "fourier":
                    phi0 = np.arccos(-float(fields[4])/(4.0*float(fields[5])))
                    kamber = float(fields[2])/2 *(-float(fields[4])*np.cos(phi0) - 4*float(fields[5])*np.cos(2*phi0))
                    phi0 = phi0 *180/np.pi
                    angles[key] = (kamber,phi0)

            # ---------------- DIHEDRALS ----------------
            elif section == "Dihedral":
                v = float(fields[1])
                sign = float(fields[2])
                periodicity = int(fields[3])

                phase = 0.0 if sign > 0 else 180.0

                atoms = map_atoms(atom_types)
                key = normalize_key(atoms)

                dihedrals[key] = (v, phase, periodicity)
                

            # ---------------- IMPROPERS ----------------
            elif section == "Improper":
                k = float(fields[1])

                atoms = map_atoms(atom_types)
                key = tuple(atoms)

                impropers[key] = k

            # ---------------- NONBONDED ----------------
            elif section == "Pair":
                epsilon = float(fields[1])
                sigma = float(fields[2])

                atom = map_atoms([atom_types[0]])[0]
                nonbond[atom] = (sigma*0.561231, epsilon)

    return {
        "masses": masses,
        "bonds": bonds,
        "angles": angles,
        "dihedrals": dihedrals,
        "impropers": impropers,
        "nonbond": nonbond,
    }


# ============================================================
# Replace frcmod entries
# ============================================================

def update_frcmod(template_file, output_file, ff):

    section = None
    output = []

    with open(template_file) as f:

        for raw in f:

            line = raw.rstrip("\n")
            stripped = line.strip()

            if stripped in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                section = stripped
                output.append(line)
                continue

            # ---------------- MASS ----------------
            if section == "MASS":

                m = re.match(r"^(\S+)", stripped)

                if m:
                    atom = m.group(1)

                    if atom in ff["masses"]:
                        mass = ff["masses"][atom]
                        output.append(f"{atom:<4} {mass:8.4f}   0.000")
                        continue

            # ---------------- BOND ----------------
            elif section == "BOND":

                m = re.match(r"^\s*([A-Za-z0-9_]+)\s*-\s*([A-Za-z0-9_]+)", stripped)

                if m:

                    a1, a2 = m.groups()
                    key = normalize_key(map_atoms([a1, a2]))

                    if key in ff["bonds"]:
                        k, req = ff["bonds"][key]

                        output.append(
                            f"{fmt_atom(a1)}-{fmt_atom(a2)} {k:8.3f} {req:8.4f}"
                        )
                        continue

            # ---------------- ANGLE ----------------
            elif section == "ANGLE":

                m = re.match(
                    r"^\s*([A-Za-z0-9_]+)\s*-\s*([A-Za-z0-9_]+)\s*-\s*([A-Za-z0-9_]+)",
                    stripped
                )

                if m:

                    a1, a2, a3 = m.groups()
                    key = normalize_key(map_atoms([a1, a2, a3]))

                    if key in ff["angles"]:
                        k, theta = ff["angles"][key]

                        output.append(
                            f"{fmt_atom(a1)}-{fmt_atom(a2)}-{fmt_atom(a3)} {k:8.3f} {theta:8.3f}"
                        )
                        continue

            # ---------------- DIHEDRAL ----------------
            elif section == "DIHE":

                m = re.match(
                    r"^\s*([A-Za-z0-9_]+)\s*-\s*([A-Za-z0-9_]+)\s*-\s*([A-Za-z0-9_]+)\s*-\s*([A-Za-z0-9_]+)",
                    stripped
                )

                if m:

                    atoms = list(m.groups())
                    key = normalize_key(map_atoms(atoms))

                    if key in ff["dihedrals"]:
                        v, phase, periodicity = ff["dihedrals"][key]

                        output.append(
                            f"{'-'.join(fmt_atom(a) for a in atoms):<20} 1 {v:8.3f} {phase:8.3f} {periodicity:4d}"
                        )
                        continue

            # ---------------- NONBON ----------------
            elif section == "NONBON":

                m = re.match(r"^(\S+)", stripped)

                if m:

                    atom = m.group(1)

                    if atom in ff["nonbond"]:
                        sigma, epsilon = ff["nonbond"][atom]

                        output.append(
                            f"{atom:<4} {sigma:10.4f} {epsilon:10.4f}"
                        )
                        continue

            output.append(line)

    with open(output_file, "w") as f:
        f.write("\n".join(output))


# ============================================================
# MAIN
# ============================================================

def main():

    if len(sys.argv) != 4:
        print("Usage: python fill_frcmod.py coeffs.txt template.frcmod output.frcmod")
        sys.exit(1)

    coeff_file, template_file, output_file = sys.argv[1:4]

    ff = parse_lammps_coeffs(coeff_file)
    update_frcmod(template_file, output_file, ff)

    print(f"Wrote updated frcmod to: {output_file}")


if __name__ == "__main__":
    main()
