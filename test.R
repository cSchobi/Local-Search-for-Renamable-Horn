#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
file = args[1]
file = "quinn"
#load .csv file
file.path("statistics", paste(file, ".csv", sep=""))
experiments = read.csv(file.path("statistics", paste(file, ".csv", sep="")), 
                       header = FALSE, sep = ",", row.names = 1)

#setup output file
png(filename=file.path("plots/", paste(file, ".png", sep="")))
#plot
col_set = rainbow(nrow(experiments))
matplot(t(experiments), type = "l", lty = "solid", col= col_set, ylab = "non Horn clauses", xlab = "flips")
legend("topright", rownames(experiments), lwd=c(2.5,2.5,2.5,2.5), col = col_set)
#close output file
dev.off()
