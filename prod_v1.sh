#!/bin/bash

# load gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

declare -i nodes
declare -f walltime 
declare -f name 

nodes=5
walltime=24:00:00


# Launch production 
for i in 6 8 ; do 
	cd n${i}_sol/production
	name=n${i}_plga 
	if [ ! -f topol.tpr ] ; then 
		cp ../../j_slurm.sh .
		c_dir="$PWD"
		sed -i "s/nodes=XX/nodes=$nodes/g" j_slurm.sh
                sed -i "s/name=XXX/name=$name/g" j_slurm.sh
                sed -i "s/-np XX/-np $((16*$nodes))/g" j_slurm.sh
                sed -i "s|workdir=XXX|workdir="$c_dir"|g" j_slurm.sh		
		sed -i "s/time=XX/time=$walltime/g" j_slurm.sh
		ln -s ../../md_prod.mdp
		ln -s ../npt_eq/confout.gro
		ln -s ../em_ions/topol.top
		ln -s ../npt_eq/restart.cpt
		ln -s ../npt_eq/md.log
		ln -s ../npt_eq/ener.edr
		ln -s /suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/bash_sub_scripts/MOD_hem_add_amber99sb-ildn.ff
		gmx_mpi grompp -f md_prod.mdp -t restart.cpt -c confout.gro -p topol.top
		rm restart.cpt
		sbatch -p ckpt -A pfaendtner-ckpt j_slurm.sh
		wait
	fi
	cd ../../
done
