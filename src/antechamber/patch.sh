#!/bin/sh

for program in am1bcc$SFX \
               antechamber$SFX \
               atomtype$SFX \
               bondtype$SFX \
               espgen$SFX \
               parmchk2$SFX \
               prepgen$SFX \
               residuegen$SFX \
               respgen$SFX; do
    cat bin_template | sed "s/REPLACE_ME/$program/" > $BINDIR/$program
    chmod +x $BINDIR/$program
done
