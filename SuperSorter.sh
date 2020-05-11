#!/bin/bash
file=$1
index=4
outfile="Out.txt"
for ((index=4; index<=8; index++)); do
	sort ${file} -k ${index} -r | cut -f 1 | head -3 | tail -2 >> ${outfile}
done
