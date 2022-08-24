#!/bin/bash

# Steps to generate .top and .gro files when setting NME and ACE cap atom charges for a molecule, after the QM Gaussian calculation

# 1) Source Amber tools
source /gscratch/pfaendtner/sarah/codes/amber20/amber20/amber.sh

# 2) Convert Gaussian ESP data into the RESP format 
espgen -o LDN.esp -i LDN.log

# 3) Generate Antechamber structure file (.ac file), must specify charge of molecule
antechamber -fi pdb -i LDN.pdb -o LDN.ac -fo ac -rn LDN -nc -2

# We are performing 2 rounds of RESP fitting:
# 	Stage 1 : All atoms, except NME or ACE caps (Set N = -1 for these atoms, order of atoms is according to pdb file), are allowed to vary.
# 	Stage 2: All Methylene or methyl hydrogens atoms are constrained to have equal charge and fit (Will depend on molecule)

# 4) Generate RESP input file for round 1
respgen -i LDN.ac -o LDN_resp1.in -f resp1 -e 1

# 5) Perform first round of RESP with the file input file to get a qout file with fit charges
resp -O -i LDN_resp1.in -o LDN_RESP1.OUT -e LDN.esp -t qout 

# 6) modify round 1 RESP input file for ACE or NME cap atoms, remove punch, esout and ##_RESP1.OUT files
vi LDN_resp1.in
# modify section under resp run #1 and add iqopt=2 (tells RESP to to read initial charges from a qin file to fix cap atom charges
# change cap atoms in second column to -1 under the 1.0 section, first column is atomic number, order of atoms is according to pdb file 
rm LDN_RESP1.OUT
rm punch 
rm esout
mv qout LDN.qin

# 7) Replace all number with zeros and replace cap atom charges with AMBER FF9X charges 
vi LDN.qin
:%s/[0-9]\+\.[0-9]\+\|[0-9]\+/0.000000/g
# Modify certain zeros (follow atom order in pdb files) for cap atom to their AMBER charges; they should sum to zero 

# 8) Run the first stage of RESP
resp -O -i LDN_resp1.in -e LDN.esp -q LDN.qin -t qout -o LDN_RESP1.OUT

# 9) Generate RESP input file for round 2
respgen -i LDN.ac -o LDN_resp2.in -f resp2 -e 1 

# 10) Modify round 2 RESP input file 
vi LDN_resp2.in
# change cap atoms in second column to -1 under the 1.0 section, first column is atomic number, order of atoms is according to pdb file

# 11) Run the second stage of RESP
resp -O -i LDN_resp2.in -e LDN.esp -q qout -t QOUT -o LDN_RESP2.OUT

# 12) Generate anterchamber file with charges from 2 rounds of RESP
antechamber -i LDN.ac -fi ac -o LDN_chrg.ac -fo ac -c rc -cf QOUT

# 13) Generate antechamber mol2 file
antechamber -fi ac -i LDN_chrg.ac -fo mol2 -o LDN.mol2 

# 14) Run parmchk 
parmchk2 -i LDN.mol2 -o LDN.frcmod -f mol2 -a Y

# 15) Run tleap to get partop and inpcrd file  
# make leap.inp file
cat << HERE >| leap.inp
source /gscratch/pfaendtner/sarah/codes/amber20/amber20/dat/leap/cmd/leaprc.GLYCAM_06j-1
source /gscratch/pfaendtner/sarah/codes/amber20/amber20/dat/leap/cmd/leaprc.gaff
loadamberparams LDN.frcmod
LDN = loadmol2 LDN.mol2
saveoff LDN LDN.lib
saveamberparm LDN LDN.partop LDN.inpcrd
quit 

HERE
# Execute tleap 
tleap -f leap.inp

# 16)  Run Acpype to get .gro and .top files 
python acpypi.py -p LDN.partop -x LDN.inpcrd
