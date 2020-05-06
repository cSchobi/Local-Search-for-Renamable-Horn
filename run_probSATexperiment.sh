#!/bin/bash

EXPERIMENTS=("probSATexp_" "probSATpoly_")
MODES=("3" "4")
CBSTARTS=("1.0" "-2.0")
CBENDS=("15.0" "10")
CMSTARTS=("0.0" "-2.0")
CMENDS=("15.0" "10.0")
SKIPS=("0" "1")
NPOINTS=20
for SKIP in ${SKIPS[@]}
do
	echo "skip=$SKIP"
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
			filename="$EXPERIMENT${filename%.*}[SKIP=$SKIP]"
			echo "storing in $filename"
			./solver $inputfile -mode $MODE -skip $SKIP > "./statistics/$filename.csv"
			./plotProbSAT.m $filename $cmStart $cmEnd $cbStart $cbEnd $NPOINTS
		done
	done
done
