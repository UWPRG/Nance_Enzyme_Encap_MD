#!/bin/bash

# load gromacs
module load icc_17-impi_2017
source /gscratch/pfaendtner/cdf6gc/codes/PRG_USE/JUNE13_2017/ARUSHI/gromacs-2016.3/bin/bin/GMXRC

counter=0

echo $PWD

#mkdir rdf_pegawmix
#echo $PWD

#cd rdf_pegawmix/

#for j in 6 8 10 20 30 ; do
#	mkdir -p n${j}peg_prodaw
#done

#cd ../
#echo $PWD



for i in 30 ; do
        cd n${i}_sol/production
        echo $PWD
        if [ -f topol.tpr ] && [ confout.gro ] ; then
		#ln -s ../../m_ndx.txt
		#ln -s ../../npeg_wat.txt
		#ln -s ../../npeg_ace.txt
		#ln -s ../../psol_rdf.txt
		#ln -s ../../OaceHwat.txt
		#ln -s ../../OcaeOcae.txt
		#ln -s ../../OwatCace.txt
		#ln -s ../../OwatOwat.txt	
		#cat m_ndx.txt | gmx_mpi make_ndx -f confout.gro 
		#cat psol_rdf.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}sol.xvg
		cat npeg_wat.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}watonly.xvg
                cat npeg_ace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o rdf_n${i}aceonly.xvg
		#cat OaceHwat.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o OaHw_n${i}aw.xvg
		#cat OaceOace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o OaOa_n${i}aw.xvg
		#cat OwatCace.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o OwCa_n${i}ace.xvg
		#cat OwatOwat.txt | gmx_mpi rdf -f traj.trr -s topol.tpr -n index.ndx -o OwOw_n${i}ace.xvg
		wait ${!}
        fi
	echo $PWD
	#cp rdf_n${i}sol.xvg  ../../rdf_pegawmix
	cp rdf_n${i}watonly.xvg ../../rdf_pegawmix
	cp rdf_n${i}aceonly.xvg ../../rdf_pegawmix
	#cp OaHw_n${i}aw.xvg ../../rdf_pegawmix
	#cp OaOa_n${i}aw.xvg ../../rdf_pegawmix
	#cp OwCa_n${i}ace.xvg ../../rdf_pegawmix
	#cp OwOw_n${i}ace.xvg ../../rdf_pegawmix
	cd ../../
        echo $PWD
        counter=$((counter+1))

done

