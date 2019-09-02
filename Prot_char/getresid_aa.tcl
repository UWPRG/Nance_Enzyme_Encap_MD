set OUTPUT [open "resid_aromatic.txt" w]
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {

	set sel [atomselect top "resname PHE TRP TYR HIS and within 4 of (resname PLG )" frame $frame]
	set resid [lsort -unique [$sel get resid]]
 	puts "$frame"
 	puts $OUTPUT $resid
}
close $OUTPUT
set OUTPUT [open "resid_hydrophobic.txt" w]
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {
 
	set sel [atomselect top "resname ALA ILE LEU VAL GLY PRO and within 4 of (resname PLG )" frame $frame]
	set resid [lsort -unique [$sel get resid]]
 	puts "$frame"
 	puts $OUTPUT $resid
}
close $OUTPUT
set OUTPUT [open "resid_polar.txt" w]
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {
 
	set sel [atomselect top "resname ASN CYS GLN MET SER THR and within 4 of (resname PLG )" frame $frame]
	set resid [lsort -unique [$sel get resid]]
 	puts "$frame"
 	puts $OUTPUT $resid
}
close $OUTPUT
set OUTPUT [open "resid_negative.txt" w]
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {
 
	set sel [atomselect top "resname ASP GLU and within 4 of (resname PLG )" frame $frame]
	set resid [lsort -unique [$sel get resid]]
 	puts "$frame"
 	puts $OUTPUT $resid
}
close $OUTPUT
set OUTPUT [open "resid_positive.txt" w]
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {
 
	set sel [atomselect top "resname ARG HIP LYS and within 4 of (resname PLG )" frame $frame]
	set resid [lsort -unique [$sel get resid]]
 	puts "$frame"
 	puts $OUTPUT $resid
}
close $OUTPUT
quit


