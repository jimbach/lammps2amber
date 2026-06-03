#!/bin/bash

[[ -f "$1.cif" ]] || { echo "Error: $1.cif not found" >&2; exit 1; }
[[ -f "$1.data" ]] || { echo "Error: $1.data not found" >&2; exit 1; }

shortname="${1:0:3}"


command -v parmchk >/dev/null 2>&1 || { echo "Error: parmchk not found in PATH" >&2; exit 1; }
command -v antechamber >/dev/null 2>&1 || { echo "Error: antechamber not found in PATH" >&2; exit 1; }
command -v obabel >/dev/null 2>&1 || { echo "Error: obabel not found in PATH" >&2; exit 1; }

obabel -iCIF $1.cif -opdb -O $1.pdb
antechamber -i $1.pdb -fi pdb -o $1.mol2 -fo mol2 -j 5 -at sybyl -dr no 

#replace the molname
sed -i "s/UNL/$shortname/g" $1.mol2

${LAMMPS2AMBERHOME}/tools/rename_mol2.py $1.data $1.mol2 renamed.mol2
mv renamed.mol2 $1.mol2


${LAMMPS2AMBERHOME}/tools/copy_bonds_to_mol2.py $1.mol2 $1.data out.mol2
mv out.mol2 $1.mol2


${LAMMPS2AMBERHOME}/tools/make_bigger_unitcell.py $1.mol2 $1.data out.mol2
mv out.mol2 $1.mol2

${LAMMPS2AMBERHOME}/tools/uniquely_label_atoms_mol2.py $1.mol2 out.mol2
mv out.mol2 $1.mol2


${LAMMPS2AMBERHOME}/tools/gen_frcmod_template.py $1.mol2 $shortname.frcmod


${LAMMPS2AMBERHOME}/tools/fill_frcmod2.py $1.data $shortname.frcmod filled.frcmod 
mv filled.frcmod $shortname.frcmod

mv removed.frcmod $shortname.frcmod

#boxx=$(grep "_cell_length_a" $1.cif | awk '{print $2}')
#boxy=$(grep "_cell_length_b" $1.cif | awk '{print $2}')
#boxz=$(grep "_cell_length_c" $1.cif | awk '{print $2}')

boxx=$(grep "xlo xhi" $1.data | awk '{print $2}')
boxy=$(grep "ylo yhi" $1.data | awk '{print $2}')
boxz=$(grep "zlo zhi" $1.data | awk '{print $2}')


filename="$1"

grep -v 0.0000 $shortname.frcmod > nozero.frcmod

echo "building a tleap file"
echo "#################################################"
tee $filename.tleap.in << END
source leaprc.uff
loadamberparams nozero.frcmod
$shortname = loadmol2 $filename.mol2
set $shortname box {$boxx $boxy $boxz}
saveoff $shortname $shortname.lib 
loadoff $shortname.lib
saveamberparm $shortname $filename.prmtop $filename.inpcrd
quit
END
echo "#################################################"
echo "running tleap"
tleap -f $filename.tleap.in | awk '!seen[$0]++' | grep " No " > MISSING-PARAMTERS.dat

${LAMMPS2AMBERHOME}/tools/add_missing_params.py MISSING-PARAMTERS.dat nozero.frcmod
mv nozero.frcmod $shortname.frcmod

echo "building a tleap file"
echo "#################################################"
tee $filename.tleap.in << END
source leaprc.uff
loadamberparams $shortname.frcmod
$shortname = loadmol2 $filename.mol2
set $shortname box {$boxx $boxy $boxz}
saveoff $shortname $shortname.lib 
loadoff $shortname.lib
saveamberparm $shortname $filename.prmtop $filename.inpcrd
quit
END
echo "#################################################"
echo "running tleap"
tleap -f $filename.tleap.in


echo "#################################################"
echo "if all went well .frcmod and .lib were generated"
echo "please revise all parameters from .frcmod and charges in .lib"
echo "then use regen-inpcrd-prmtop to generate the proper input for cp2k"


echo "box: "$boxx" "$boxy" "$boxz



echo "#################################################"
echo "############# autotopology finished #############"
echo "#################################################"
