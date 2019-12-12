#!/bin/bash

# load gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

counter=0

echo $PWD

mkdir plga_awmix
cd plga_awmix/
echo $PWD


for j in 6 8 10 20 30 ; do
	mkdir -p n${j}plga_awmix
done

cd ../
echo $PWD

# Launch production 
for i in 6 8 10 20 30 ; do
        cd n${i}_sol/production
        echo $PWD
#        if [ -f topol.tpr ] && [ confout.gro ] ; then
#		ln -s ../../st_s.txt
#		cat st_s.txt | gmx_mpi trjconv -s topol.tpr -f traj.trr -o noawmix_n${i}plga.xtc -pbc mol -ur compact
#		cat st_s.txt | gmx_mpi trjconv -s topol.tpr -f traj.trr -o n${i}plgaonly_awmix.pdb -dump 0	
#                wait ${!}
#        fi
#	echo $PWD	
	cp n${i}plgaonly_awmix.pdb ../../plga_awmix/n${i}plga_awmix
	cp ener.edr ../../plga_awmix/n${i}plga_awmix
	cp md.log ../../plga_awmix/n${i}plga_awmix
	cp topol.top ../../plga_awmix/n${i}plga_awmix
	cp topol.tpr ../../plga_awmix/n${i}plga_awmix
	cp confout.gro  ../../plga_awmix/n${i}plga_awmix	
	cp noawmix_n${i}plga.xtc ../../plga_awmix/n${i}plga_awmix	
	cd ../../
        echo $PWD
        counter=$((counter+1))

done

