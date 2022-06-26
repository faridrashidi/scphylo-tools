# BiocManager::install("copynumber")
# install.packages("sequenza")

args <- commandArgs(trailingOnly=TRUE)
library(sequenza)
out_dir <- paste0(args[1],"/",args[2],".sequenza")

test <- sequenza.extract(
    paste0(args[1],"/",args[2],".50.seqz.gz"),
    chromosome.list=paste("chr",c((1:22),"X","Y"),sep=""),
    verbose=FALSE
)
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract=test, cp.table=CP, sample.id=args[2], out.dir=out_dir)

cint <- get.ci(CP)
cellularity <- cint$max.cellularity
write.table(cellularity, paste(out_dir,paste(args[2],"_cellularity.txt",sep=""),sep="/"), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
ploidy <- cint$max.ploidy
write.table(ploidy, paste(out_dir,paste(args[2],"_ploidy.txt",sep=""),sep="/"), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

sequenza2PyClone <- function(mut.tab, seg.cn, sample.id, norm.cn=2) {
    mut.tab <- cbind(mut.tab[,c("chromosome","position","good.reads","F","mutation")], CNt=NA, A=NA, B=NA)
    for (i in 1:nrow(seg.cn)) {
       pos.filt <- mut.tab$chromosome == seg.cn$chromosome[i] & mut.tab$position >= seg.cn$start.pos[i] & mut.tab$position <= seg.cn$end.pos[i]
       mut.tab[pos.filt, c("CNt", "A", "B")] <- seg.cn[i, c("CNt", "A", "B")]
    }
    id          <- paste(sample.id, mut.tab$chromosome, mut.tab$position, sep=":")
    var.counts  <- round(mut.tab$good.reads * mut.tab$F, 0)
    nor.counts  <- mut.tab$good.reads - var.counts
    pyclone.tsv <- data.frame(
        mutation_id=id, ref_counts=nor.counts, var_counts=var.counts,
        normal_cn=norm.cn, minor_cn=mut.tab$B, major_cn=mut.tab$A,
        variant_case=sample.id, variant_freq=mut.tab$F, genotype=mut.tab$mutation
    )
    pyclone.tsv <- pyclone.tsv[pyclone.tsv$major_cn != 0, ]
    na.exclude(pyclone.tsv)
}
mut.tab <- read.table(paste(out_dir,paste(args[2],"_mutations.txt",sep=""),sep="/"), header=TRUE)
seg.cn <- read.table(paste(out_dir,paste(args[2],"_segments.txt",sep=""),sep="/"), header=TRUE)
pyclone.tsv <- sequenza2PyClone(
    mut.tab=mut.tab,
    seg.cn=seg.cn,
    sample.id=args[2]
)
write.table(pyclone.tsv, paste(out_dir,paste(args[2],"_pyclone.tsv",sep=""),sep="/"), row.names=FALSE, quote=FALSE, sep="\t")
