#!/bin/bash

for file in examples/*.cnf
do
	echo "processing $file"
	filename=$(basename $file)
	filename="boxPlot${filename%.*}"
	echo "storing in $filename"
	./solver $file -mode 5 > "./statistics/$filename.csv"
	Rscript plotBoxPlot.R $filename
done
