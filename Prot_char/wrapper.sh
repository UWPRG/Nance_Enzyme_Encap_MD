#!/bin/bash

# Remove files previously generated from running this bash file 
rm -fr resid.txt resid2.txt occ.txt

#Execute the getresi.tcl script
vmd -dispdev text ../confout.gro ../plga_25ns_p1.xtc < getresid.tcl

# for each line in resid.txt, organize each residue number in descending order and output resid2.txt 
sed 's/ /\n/g' resid.txt | awk 'NF > 0' | sort -nk 1 > resid2.txt

count=`wc -l resid.txt | awk '{print $1}'`
for ((i=1;i<=583;i++)); do
 	echo $count	
	awk -v r=$i -v count=$count 'BEGIN{occ=0}{if($1==r){occ++}}END{print r,occ}' resid2.txt >> occ_25ns_p1.txt
done;
