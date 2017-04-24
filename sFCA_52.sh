#!/bin/bash

date=20161113
begin=10001  #steeds plus 8!
counter=0

MEAN_Z=9.4 #um from volume

np=1
MAXPROCESSOR=12

total_line=200 

line=1
#run through cell areas file
while [ $line -le $total_line ] ; do
    #get cell Area from file
    Area=`sed -n $line'p' ColFRI_Potsdam_uns_linux_20161006.txt`
  
    line=$(($line+1))
    Vol=`echo "scale=4; $MEAN_Z*(10^(-3))*$Area" | tr -d $'\r' | bc` 
    #tr gets rid of weird unwritten character
  
   
    for x in {1..50}
    do 
	VAR_file=$(($begin + $counter))
	time ./FCA_52 $date $VAR_file $Vol > $date"_z"$VAR_file.txt &
	
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










