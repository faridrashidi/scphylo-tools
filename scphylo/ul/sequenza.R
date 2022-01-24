# BiocManager::install("copynumber")
# install.packages("sequenza")

args <- commandArgs(trailingOnly=TRUE)
library(sequenza)
test <- sequenza.extract(paste0(args[1],"/",args[2],".50.seqz.gz"), verbose=TRUE)
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract=test, cp.table=CP, sample.id=args[2], out.dir=paste0(args[1],"/",args[2],".sequenza"))
