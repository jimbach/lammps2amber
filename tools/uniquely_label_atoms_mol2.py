#!/usr/bin/env python3

import re
import sys
from collections import defaultdict


#def element(name: str) -> str:
#    """Extract element prefix (C from C15, Zn from Zn2, etc.)."""
#    m = re.match(r"([A-Za-z]+)", name.strip())
#    return m.group(1) if m else name.strip()

def element(name: str) -> str:
    """
    Extract element symbol.

    Zn → Z
    Fe → F
    Cl → C
    C15 → C
    H2  → H
    """

    m = re.match(r"([A-Za-z]+)", name.strip())
    if not m:
        return name.strip()

    el = m.group(1)

    # collapse multi-letter elements to first letter only
    if len(el) > 1:
        return el[0].upper()

    return el.upper()

def make_label(el: str, i: int) -> str:
    """
    Naming scheme:

    0          -> H
    1–9        -> H1–H9
    10–243     -> H1A–H9Z   (234)
    244–477    -> HA1–HZ9   (234)
    478–1153   -> HAA–HZZ   (676)
    """

    if i == 0:
        return el

    if i <= 9:
        return f"{el}{i}"

    i -= 10  # now start extended blocks

    # Block 1: H1A–H9Z (9 × 26)
    block1 = 9 * 26  # 234
    if i < block1:
        d = i // 26 + 1
        l = i % 26
        return f"{el}{d}{chr(65 + l)}"

    i -= block1

    # Block 2: HA1–HZ9 (26 × 9)
    block2 = 26 * 9  # 234
    if i < block2:
        l = i // 9
        d = i % 9 + 1
        return f"{el}{chr(65 + l)}{d}"

    i -= block2

    # Block 3: HAA–HZZ (26 × 26)
    l1 = i // 26
    l2 = i % 26

    if l1 >= 26:
        raise ValueError(f"Too many atoms for element {el} (max 1154)")

    return f"{el}{chr(65 + l1)}{chr(65 + l2)}"


def relabel_mol2(inp, outp):

    counts = defaultdict(int)
    in_atom = False

    with open(inp, "r") as fin, open(outp, "w") as fout:

        for line in fin:

            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                counts.clear()
                fout.write(line)
                continue

            if in_atom and line.startswith("@<TRIPOS>"):
                in_atom = False
                fout.write(line)
                continue

            if not in_atom:
                fout.write(line)
                continue

            if not line.strip():
                fout.write(line)
                continue

            # Fixed-width MOL2 assumption:
            # 0–6   atom id (7 chars)
            # 7–13  atom name (7 chars)
            # 14+   rest
            atom_id = line[:7]
            old_name = line[7:14]
            rest = line[14:]

            el = element(old_name)
            idx = counts[el]
            counts[el] += 1

            new_name = make_label(el, idx)

            if len(new_name) > 7:
                raise ValueError(
                    f"Atom name '{new_name}' exceeds 7-character field"
                )

            fout.write(atom_id + f"{new_name:>7}" + rest)


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.mol2 output.mol2")
        sys.exit(1)

    relabel_mol2(sys.argv[1], sys.argv[2])
