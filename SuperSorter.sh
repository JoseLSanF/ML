#!/bin/bash
file=$1
index=4
outfile="ChosenProteins.txt"
for ((index=4; index<=8; index++)); do
	name=$(cat ${file} | head -1 | cut -f ${index}
	echo ${name} >> ${outfile}
	sort ${file} -k ${index} -r | cut -f 1 | head -3 | tail -2 >> ${outfile}
done
