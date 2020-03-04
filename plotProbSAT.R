#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
file = args[1]
nPoints = 30
startValue = 1.0
endValue = 20.0
#load .csv file
file.path("statistics", paste(file, ".csv", sep=""))
experiments = read.csv(file.path("statistics", paste(file, ".csv", sep="")), 
                       header = FALSE, sep = ",", row.names = 1)
# setup x values
x = matrix(seq(startValue, endValue, length.out = nPoints), ncol = 1)
#setup output file
png(filename=file.path("plots/", paste(file, ".png", sep="")))
#plot
col_set = rainbow(nrow(experiments))
matplot(x, t(experiments), type = "l", lty = "solid", col= col_set, ylab = "non Horn clauses", xlab = "cb")
legend("right", rownames(experiments), lwd=c(2.5), col = col_set)
#close output file
dev.off()
