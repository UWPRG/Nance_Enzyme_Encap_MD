#!/bin/bash

# load gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

#make directories for polymer/acetone or polymer/DMSO only simulations
#counter=0

declare -a bsize
bsize=(5.5 5.5 6.5 8.0 10.0)

declare -a n_size
n_size=(5.7 5.7 6.7 8.2 10.2)

declare -a mol_n
mol_n=(1354 1354 2235 4166 8138)

for i in 30 ; do	
	echo $PWD
#	if [ ! -d n${i}_sol ] ; then 
#		mkdir -p n${i}_sol
#		mkdir -p n${i}_sol/em_ions
#		mkdir -p n${i}_sol/production
#		mkdir -p n${i}_sol/nvt_eq
#		mkdir -p n${i}_sol/npt_eq	
#		mkdir -p n${i}_sol/scaling
#		cp minim.mdp n${i}_sol/em_ions
#		cp nvt.mdp n${i}_sol/nvt_eq
#		cp npt.mdp n${i}_sol/npt_eq
#		cp md_prod.mdp n${i}_sol/production	
#		cp n${i}_plga.pdb n${i}_sol/em_ions
#		cp acetone.pdb n${i}_sol/em_ions
#		cp packmol.inp n${i}_sol/em_ions	
#	fi

#generate gro and top files, 
#before running this step, make sure the polymer pdb is in the right folders, ions don't need to be added 
	
	#launch energy minization
#	if [ ! -f n${i}_sol/em_ions/box_n${i}.gro ] && [ ! -f n${i}_sol/em_ions/posre.itp ] ; then 
#		cd n${i}_sol/em_ions
#		echo $PWD
#		cp ../../packmol.inp .
#		bcen=$(bc <<< "scale=1; ${bsize[counter]}*5")
#		btot=$(bc <<< "scale=1; ${bsize[counter]}*10") 
#		cat > packmol.inp << EOF1
##
## polymer and acetone
#
#tolerance 2.0
#filetype pdb
#output n${i}_processed.pdb
#
#structure n${i}_plga.pdb
#  number 1
#  center
#  fixed $bcen $bcen $bcen $bcen $bcen $bcen
#end structure
#
#structure acetone.pdb
#  number ${mol_n[counter]}
#  inside box 0. 0. 0. $btot $btot $btot
#end structure
#EOF1
#		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/MOD_hem_add_amber99sb-ildn.ff
#		ln -s ../../../pdb2gmx.txt 
#		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/polymer_force_field/conf_data/residuetypes.dat
#		packmol < packmol.inp 
#		gmx_mpi editconf -f n${i}_processed.pdb -box ${n_size[counter]} ${n_size[counter]} ${n_size[counter]} -o box_n${i}.gro
#		cat pdb2gmx.txt | gmx_mpi pdb2gmx -f n${i}_plga.pdb  
#		echo "acetone           ${mol_n[counter]}" >> topol.top
#		sed -i '23i\; Include acetone parameters' topol.top 
#		sed -i '24i\#include "/suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/Acetone_H2O/acetone.itp"' topol.top
#		gmx_mpi grompp -f minim.mdp -c box_n${i}.gro -p topol.top	
#		gmx_mpi mdrun -v -deffnm topol &> log.txt
#		wait ${!}
#	fi 
#
#	cd ../ 
#	echo $PWD
#generate gro and top files, 
#before running this step, make sure the polymer pdb is in the right folders, ions don't need to be added 
	cd n${i}_sol/
	echo $PWD
	#launch NVT equilibration 
	if [ ! -f nvt_eq/confout.gro ] ; then
		cd nvt_eq/
		echo $PWD
		ln -s ../em_ions/topol.top
		ln -s ../em_ions/topol.gro
		ln -s ../em_ions/posre.itp
		ln -s ../em_ions/topol.edr
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/MOD_hem_add_amber99sb-ildn.ff
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/polymer_force_field/conf_data/residuetypes.dat 
		gmx_mpi grompp -f nvt.mdp -c topol.gro -p topol.top
		gmx_mpi mdrun -cpi restart -cpo restart -append -cpt 1.0 &> log.txt
		wait ${!}
	fi

	cd ../
	echo $PWD

	#launch NPT equilibration 
	if [ ! -f npt_eq/confout.gro ] ; then
		cd npt_eq/
		echo $PWD
		ln -s ../nvt_eq/confout.gro
		ln -s ../nvt_eq/topol.top
		ln -s ../nvt_eq/restart.cpt
		ln -s ../nvt_eq/ener.edr
		ln -s ../nvt_eq/md.log
		ln -s ../em_ions/posre.itp
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/polymer_force_field/conf_data/residuetypes.dat
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/MOD_hem_add_amber99sb-ildn.ff
		gmx_mpi grompp -f npt.mdp -c confout.gro -t restart.cpt -p topol.top
		rm restart.cpt
		gmx_mpi mdrun -cpi restart -cpo restart -append -cpt 1.0 &> log.txt
		wait ${!}
	fi

	cd ../
	echo $PWD
	# Perform scaling
	cd scaling/
	echo $PWD
	
	for j in 1 2 4 5 6 ; do
		mkdir scal_n${j}
		cd scal_n${j}
		cp ../../../j_slurm.sh . 
		echo $PWD
		c_dir="$PWD"
		sed -i "s/nodes=XX/nodes=${j}/g" j_slurm.sh
		sed -i "s/name=XXX/name=scal_${j}/g" j_slurm.sh
		sed -i "s/-np XX/-np $((16*${j}))/g" j_slurm.sh
		sed -i "s|chdir=XXX|chdir="$c_dir"|g" j_slurm.sh
		ln -s ../../../md_scaling.mdp
		ln -s ../../npt_eq/confout.gro
		ln -s ../../em_ions/topol.top
		ln -s ../../npt_eq/restart.cpt
		ln -s ../../npt_eq/md.log
		ln -s ../../npt_eq/ener.edr
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/polymer_force_field/conf_data/residuetypes.dat
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/MOD_hem_add_amber99sb-ildn.ff
		gmx_mpi grompp -f md_scaling.mdp -t restart.cpt -c confout.gro -p topol.top
		rm restart.cpt
		sbatch -p ckpt -A pfaendtner-ckpt j_slurm.sh
		wait
		cd ../
		echo $PWD
	done

	cd ../../
	echo $PWD
	counter=$((counter+1))
done 

