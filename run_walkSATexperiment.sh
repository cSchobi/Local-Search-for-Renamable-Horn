#!/bin/bash

SKIPS=("0" "1")
for SKIP in ${SKIPS[@]}
do
	echo "Skip=$SKIP"
	FILE="walkSATexperiment[SKIP=$SKIP]"
	rm "./statistics/$FILE.csv"
	for file in ./examples/*.cnf
	do
		echo "processing $file"
		filename=$(basename $file)
		printf "${filename%.*}" >> "./statistics/$FILE.csv"
		./solver $file -mode 2 -skip $SKIP >> "./statistics/$FILE.csv"
	done
	Rscript plotWalkSAT.R $FILE

	FILE="WalkSATWithMakeexperiment[SKIP=$SKIP]"
	rm "./statistics/$FILE.csv"

	for file in ./examples/*.cnf
	do
		echo "processing $file"
		filename=$(basename $file)
		printf "${filename%.*}" >> "./statistics/$FILE.csv"
		./solver $file -mode 6 -skip $SKIP >> "./statistics/$FILE.csv"
	done
	Rscript plotWalkSAT.R $FILE
done
