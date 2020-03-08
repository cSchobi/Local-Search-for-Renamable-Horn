#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
fileName = args[1]
nPoints = 20
#load .csv file
experiments = read.csv(file.path("statistics", paste(fileName, ".csv", sep="")), 
                       header = FALSE, sep = ",", row.names = 1)
# setup x values
x = matrix(seq(0, 1, length.out = nPoints), ncol = 1)
#setup output file
png(filename=file.path("plots/", paste(fileName, ".png", sep="")))
#plot
col_set = rainbow(nrow(experiments))
matplot(x, t(experiments), type = "l", lty = "solid", col= col_set, ylab = "non Horn clauses", xlab = "probability of random walk")
legend("right", rownames(experiments), lwd=c(2.5), col = col_set)
#close output file
dev.off()
