#!/bin/bash

for file in ./examples/*.cnf
do
	echo "processing $file"
	filename=$(basename $file)
	./solver $file > "./statistics/${filename%.*}.csv"
	Rscript test.R ${filename%.*}
done
