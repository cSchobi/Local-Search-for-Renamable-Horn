#!/bin/bash

for file in ./examples/*.cnf
do
	echo "processing $file"
	filename=$(basename $file)
	./solver $file > "./statistics/${filename%.*}.csv"
	Rscript plotNHornTrend.R ${filename%.*}
done
