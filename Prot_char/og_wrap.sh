#!/bin/bash

rm -fr resid.txt resid2.txt occ.txt
vmd -dispdev text ../confout.gro ../pla_last100ns.xtc < getresid.tcl
sed 's/ /\n/g' resid.txt | awk 'NF > 0' | sort -nk 1 > resid2.txt

count=`wc -l resid.txt | awk '{print $1}'`
for ((i=1;i<=583;i++)); do
 	echo $count	
	awk -v r=$i -v count=$count 'BEGIN{occ=0}{if($1==r){occ++}}END{print r,occ/count}' resid2.txt >> occ.txt
done;
