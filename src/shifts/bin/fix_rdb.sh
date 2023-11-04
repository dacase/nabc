#!/bin/sh

#  fix up the rdb file made by shifts: remove spaces, sort on
#  res,atomname;  over-write input file

head -2 $1 > $1.tmp
awk 'BEGIN{OFS="\t"}NR>2{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $1 | sort -k1n -k3 >> $1.tmp
/bin/mv $1.tmp $1

