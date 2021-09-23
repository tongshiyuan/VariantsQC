#!/bin/bash

infile=$1
ref=$2
outfile=${infile}.seq
cat ${infile} | while read  A B C;
do
    echo -e -n "$A\t$B\t$C\t" && samtools faidx $ref "chr${A}:${B}-${C}" | grep -v '^>'

done >> ${outfile}