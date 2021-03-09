#!/bin/bash

# Test OL3 Chi parameters for RNA.

TESTtleap="../../../bin/tleap"
TESTsander="../../../bin/msander"
DACDIF="../../dacdif"

printf "\nTest OL3 Chi parameters for RNA.\n\n"
cat > leap.in <<EOF
source leaprc.RNA.OL3
source leaprc.water.tip3p
# TIP3P ions
loadamberparams frcmod.ionsjc_tip3p
m = loadpdb rGACC.pdb
solvateoct m TIP3PBOX 15.54988860294556361621 0.9
addions m Na+ 1
addions m Na+ 1
addions m Na+ 1
saveamberparm m rGACC.tip3p.parm7 rGACC.nomin.rst7
quit
EOF

$TESTtleap -f leap.in > leap.out

$DACDIF rGACC.tip3p.parm7.save rGACC.tip3p.parm7
$DACDIF rGACC.nomin.rst7.save rGACC.nomin.rst7

/bin/rm -f  leap.in leap.out leap.log

exit 0
