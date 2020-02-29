#!/bin/bash


FILE=walkSATexperiment
rm "./statistics/$FILE.csv"
for file in ./examples/*.cnf
do
	echo "processing $file"
	filename=$(basename $file)
	printf "${filename%.*}" >> "./statistics/$FILE.csv"
	./solver $file 2 >> "./statistics/$FILE.csv"
done
Rscript plotWalkSAT.R $FILE
