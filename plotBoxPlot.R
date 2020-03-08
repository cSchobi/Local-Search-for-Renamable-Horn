#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
fileName = args[1]
experiments = read.csv(file.path("statistics", paste(fileName, ".csv", sep="")), 
                         header = FALSE, sep = ",", row.names = 1)
png(filename=file.path("plots/", paste(fileName, ".png", sep="")))
par(mar=c(10,4,1,1))
boxplot(t(experiments), las=2, ylab = "Number of non horn clauses")
dev.off()

