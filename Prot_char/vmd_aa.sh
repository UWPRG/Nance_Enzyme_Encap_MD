#!/bin/bash

rm -fr resid_*.txt

vmd -dispdev text ../confout.gro ../plga_25ns_p1.xtc < getresid_aa.tcl

# for each line in resid.txt, organize each residue number in descending order and output resid2.txt 
#sed 's/ /\n/g' resid_aromatic.txt | awk 'NF > 0' | sort -nk 1 > resid2aro.txt

#sed 's/ /\n/g' resid_hydrophobic.txt | awk 'NF > 0' | sort -nk 1 > resid2hydro.txt

#sed 's/ /\n/g' resid_polar.txt | awk 'NF > 0' | sort -nk 1 > resid2polr.txt

#sed 's/ /\n/g' resid_negative.txt | awk 'NF > 0' | sort -nk 1 > resid2neg.txt

#sed 's/ /\n/g' resid_positive.txt | awk 'NF > 0' | sort -nk 1 > resid2pos.txt



