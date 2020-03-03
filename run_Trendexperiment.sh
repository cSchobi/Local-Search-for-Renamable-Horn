#!/bin/bash

SKIP=("0" "1")
for file in ./examples/*.cnf
do
	echo "processing $file"
	for SKIP_MODE in ${SKIP[@]}
	do
		echo "skip mode $SKIP_MODE"
		filename=$(basename $file)
		./solver $file -mode 1 -skip $SKIP_MODE > "./statistics/${filename%.*}[SKIP=$SKIP_MODE].csv"
		Rscript plotTrend.R "${filename%.*}[SKIP=$SKIP_MODE]"
	done
done

#show all runs in one plot for different skip values
#for file in ./examples/*.cnf
#do
#	echo "processing $file"
#		echo "skip mode $SKIP_MODE"
#		filename=$(basename $file)
#		./solver $file -mode 5 > "./statistics/${filename%.*}.csv"
#		Rscript plotTrend.R "${filename%.*}"
#done
