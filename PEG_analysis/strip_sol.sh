#!/bin/bash

# load gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

counter=0

echo $PWD

mkdir tout_pegwater
cd tout_pegwater/
echo $PWD

for j in 6 8 10 20 30 ; do
	mkdir -p n${j}peg_prod
done

cd ../
echo $PWD



for i in 6 8 10 20 30 ; do
        cd n${i}_sol/production
        echo $PWD
#        if [ -f topol.tpr ] && [ confout.gro ] ; then
#		ln -s ../../st_s.txt
#		cat st_s.txt | gmx_mpi trjconv -s topol.tpr -f traj.trr -o nowat_n${i}peg.xtc -pbc mol -ur compact
#		cat st_s.txt | gmx_mpi trjconv -s topol.tpr -f traj.trr -o n${i}pegonly.pdb -dump 0	
#                wait ${!}
#        fi
#	echo $PWD	
#	cp n${i}pegonly.pdb ../../tout_pegwater/n${i}peg_prod
#	cp ener.edr ../../tout_pegwater/n${i}peg_prod
#	cp md.log ../../tout_pegwater/n${i}peg_prod
#	cp topol.top ../../tout_pegwater/n${i}peg_prod
#	cp topol.tpr ../../tout_pegwater/n${i}peg_prod
	cp confout.gro  ../../tout_pegwater/n${i}peg_prod	
	cd ../../
        echo $PWD
        counter=$((counter+1))

done

