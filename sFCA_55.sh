#!/bin/bash

date=20170222
begin=440001  #steeds plus 8!
counter=0

N_loci=4
MEAN_Z=9.4 #um from volume

np=1
MAXPROCESSOR=12

total_line=200 
for y in $(seq 1 $N_loci)
do 
    line=1
    #run through cell areas file
    while [ $line -le $total_line ] ; do
        #get cell Area from file
        #Area=`sed -n $line'p' ColFRI_mRNA_uns_linux_20160929.txt`
	Area=`sed -n $line'p' ColFRI_Potsdam_uns_linux_20161006.txt`
        #echo $Area
	line=$(($line+1))
	Vol=`echo "scale=4; $MEAN_Z*(10^(-3))*$Area" | tr -d $'\r' | bc` 
        #tr gets rid of weird unwritten character
	
	for x in {1..50}
	do 
	    VAR_file=$(($begin + $counter))
	    time ./FCA_55 $date $VAR_file $Vol $N_loci > $date"_z"$VAR_file.txt & #extra command: nohup
	
	    ((counter += 1))
	    if test $np -ge $MAXPROCESSOR
	    then
		wait
		np=1
	    else
		let np=$np+1
	    fi
	done
    done
done

: << EOF

EOF










