# Install DESeq2 and data
#renv::init()
#renv::install("bioc::DESeq2")
#renv::install("bioc:apeglm") # For 
#renv::install(c("tidyverse", "bioc::tximport", "bioc::tximportData"))
#renv::install("pheatmap")
#renv::install("bioc::vsn", "hexbin)

library("tximport")
library("readr")
library("tximportData")
# Import data
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir, "samples.txt"), header=TRUE)
samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]

# Prepare file to load
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

# import
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# Use
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Prefiltering
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep,]

dds <- DESeq(ddsTxi)

res <- results(dds, contrast=c("condition", "B", "A"))
resLFC <- lfcShrink(dds, coef="condition_B_vs_A", type="apeglm")
plotCounts(dds, gene="ENSG00000000460.16", intgroup="condition")

plotMA(res)
plotMA(resLFC)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(vsd)["condition"])
pheatmap(as.matrix(assay(vsd)[select,]),  cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE,annotation_col = df)


# PCA
plotPCA(vsd, intgroup=c("condition"))

# More info in https://github.com/hbctraining/DGE_workshop/tree/master/lessons
