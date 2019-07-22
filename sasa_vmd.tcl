set OUTPUT [open "prot_chA.txt" w]

# Get list of residues 
set allsel [atomselect top "protein chain A"]
set residlist [lsort -unique [$allsel get resid]]

# Make atom selection, calculate sasa for each residue 

for each r $residlist {
	set sel [atomselect top "resid $r"]
	set rsasa [measure sasa 1.4 $allsel -restrict $sel]
	puts $OUTPUT "residue: $r, sasa: $rsasa"
}
close $OUTPUT
quit
