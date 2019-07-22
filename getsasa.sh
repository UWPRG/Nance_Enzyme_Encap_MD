#!/bin/bash

# Calculate Solvent Accessible Surface Area
# Source: https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/4190.html
# Vmd code is in sasa_vmd.tcl

#export PATH=$PATH:$HOME/vmd_1.9.4/VMD\ 1.9.4.app/Contents/vmd/vmd_MACOSXX86
#vmd

$HOME/vmd_1.9.4/VMD\ 1.9.4.app/Contents/vmd/vmd_MACOSXX86 -dispdev text cat_btprot.pdb < sasa_vmd.tcl

