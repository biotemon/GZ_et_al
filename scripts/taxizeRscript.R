#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("/Users/tito-admin/Tito/JOYELABACKUP/GZ_BACKUP/Temp14/get_taxonomy/")
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

library("taxize")

filename <- args[1]
dat <- read.csv(file = filename, header = FALSE)
vec <- as.vector(dat[,1])

out <- classification(vec, db = 'itis')
#out <- classification(vec, db = 'ncbi')
print(out[])
