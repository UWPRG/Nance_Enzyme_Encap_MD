set OUTPUT [open "resid.txt" w]
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {
 
	set sel [atomselect top "protein and within 4 of ( resname PLG )" frame $frame]
	set resid [lsort -unique [$sel get resid]]
 	puts "$frame"
 	puts $OUTPUT $resid
}
close $OUTPUT
quit
