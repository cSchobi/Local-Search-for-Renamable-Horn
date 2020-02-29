#!/bin/bash

for file in ./examples/*.cnf
do
	echo "processing $file"
	filename=$(basename $file)
	./solver $file 1 > "./statistics/${filename%.*}.csv"
	Rscript plotTrend.R ${filename%.*}
done
