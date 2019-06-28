
counter=0

declare -a bsize
bsize=(5.5 5.5 6.5 8.0 10.0)

declare -a n_size
n_size=(1292 1292 2131 3975 7764)

for i in 6 8 ; do
        echo ${bsize[counter]}
        echo ${i}
	bcen=$(bc <<< "scale=1; ${bsize[counter]}*5")
	echo $bcen
	btot=$(bc <<< "scale=1; ${bsize[counter]}*10")
	echo $btot
	cat > packmol.inp << EOF1
#
# polymer and acetone

tolerance 2.0
filetype pdb
output n${i}_processed.pdb

structure n${i}_plga.pdb
  number 1
  center
  fixed $bcen $bcen $bcen $bcen $bcen $bcen
end structure

structure acetone.pdb
  number ${n_size[counter]}
  inside box 0. 0. 0. $btot $btot $btot
end structure
EOF1
	counter=$((counter+1))
	echo $counter
done
