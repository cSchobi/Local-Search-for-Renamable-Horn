#!/bin/bash

EXPERIMENTS=("probSATexp_" "probSATpoly_")
MODES=("3" "4")
CBSTARTS=("1.0" "1.0")
CBENDS=("4.0" "4.0")
CMSTARTS=("0.0" "-1.0")
CMENDS=("2.0" "1.0")
NPOINTS=20
for ((i=0; i < ${#EXPERIMENTS[@]};++i));do
	EXPERIMENT=${EXPERIMENTS[i]}
	MODE=${MODES[i]}
	cmStart=${CMSTARTS[i]}
	cmEnd=${CMENDS[i]}
	cbStart=${CBSTARTS[i]}
	cbEnd=${CBENDS[i]}
	echo "processing $EXPERIMENT"
	echo "mode $MODE"
	for inputfile in ./examples/*.cnf
	do
		echo "processing $inputfile"
		filename=$(basename $inputfile)
		filename=$EXPERIMENT${filename%.*}
		echo "storing in $filename"
		./solver $inputfile -mode $MODE > "./statistics/$filename.csv"
		./plotProbSAT.m $filename $cmStart $cmEnd $cbStart $cbEnd $NPOINTS
	done
done
