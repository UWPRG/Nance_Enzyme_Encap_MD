#!/bin/bash

# load gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

counter=0

echo $PWD

#mkdir rdfplga_awmix
#echo $PWD

#for j in 6 8 10 20 30 ; do
#	mkdir -p n${j}peg_prod
#done

#cd ../
#echo $PWD



for i in 30 ; do
        cd n${i}_sol/production
        echo $PWD
        if [ -f topol.tpr ] && [ confout.gro ] ; then
		#ln -s ../../m_ndx.txt
		#ln -s ../../OO_dmso.txt
		#ln -s ../../pbb_sol.txt
		#ln -s ../../pbb_water.txt
		#ln -s ../../pbb_ace.txt
		#ln -s ../../pbb_sdmso.txt
		#ln -s ../../pcbb_oace.txt
		#ln -s ../../pcbb_Cace.txt
		#ln -s ../../pobb_owat.txt
		#ln -s ../../pcbb_owat.txt
		#ln -s ../../pobb_Cace.txt
		#cat m_ndx.txt | gmx_mpi make_ndx -f confout.gro 
		#cat pbb_water.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}watonly.xvg
                #cat pbb_ace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}aceonly.xvg
		#cat pbb_sol.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}sol.xvg
		#cat pcbb_oace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}Cbb_O1ace.xvg
		#cat pcbb_Cace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}Cbb_C1ace.xvg
		#cat pobb_owat.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}Obb_OWwat.xvg
		#cat pcbb_owat.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}Cbb_OWwat.xvg
		cat pobb_Cace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}Obb_C1ace.xvg
		#cat OO_dmso.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o OO_n${i}dmso.xvg	
		wait ${!}
        fi
	echo $PWD
	cp rdf_n${i}watonly.xvg  ../../rdfplga_awmix
	cp rdf_n${i}aceonly.xvg ../../rdfplga_awmix
	cp rdf_n${i}Cbb_O1ace.xvg ../../rdfplga_awmix
	cp rdf_n${i}sol.xvg  ../../rdfplga_awmix
	cp rdf_n${i}Cbb_C1ace.xvg  ../../rdfplga_awmix
	cp rdf_n${i}Obb_OWwat.xvg  ../../rdfplga_awmix
	cp rdf_n${i}Cbb_OWwat.xvg  ../../rdfplga_awmix
	cp rdf_n${i}Obb_C1ace.xvg  ../../rdfplga_awmix
	#cp OO_n${i}dmso.xvg ../../rdfplga_awmix
	cd ../../
        echo $PWD
        #counter=$((counter+1))

done

