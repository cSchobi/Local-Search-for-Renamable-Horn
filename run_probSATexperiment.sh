#!/bin/bash

FILES=("probSATexp_experiment" "probSATpoly_experiment")
MODES=("3" "4")
for ((i=0; i < ${#FILES[@]};++i));do
	FILE=${FILES[i]}
	MODE=${MODES[i]}
	echo "processing $FILE"
	echo "mode $MODE"
	rm "./statistics/$FILE.csv"

	for inputfile in ./examples/*.cnf
	do
		echo "processing $inputfile"
		inputfilename=$(basename $inputfile)
		printf "${inputfilename%.*}" >> "./statistics/$FILE.csv"
		./solver $inputfile -mode $MODE >> "./statistics/$FILE.csv"
	done
	Rscript plotProbSAT.R $FILE
done
