# Check if ggplot2 is installed, if not installed create a personal library and install there, first time installation may take time!.
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
packages <- c("ggplot2")
install.packages(setdiff(packages, rownames(installed.packages())))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE))
	BiocManager::install("DESeq2")

library(ggplot2)
#suppressMessages(library(DESeq2, quietly = T))
suppressMessages(library(DESeq2))
#suppressPackageStartupMessages(library(DESeq2, quietly = T))
options(warn=-1)
args <- commandArgs(TRUE)
countData <- read.csv(file=args[1], header = TRUE, sep = ",")
metaData <- read.csv(file=args[2], header = TRUE, sep = ",")


dds <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design=~group, tidy = TRUE)
dds <- DESeq(dds,  fitType='local', quiet=TRUE)
res <- results(dds)
res <- res[order(res$padj),]

pdf(file=args[3], width=8,height=8)
#pdf("volcano_plot.pdf", width=8,height=8)
par(mfrow=c(2,3))
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
garbage <- dev.off()

vsdata <- varianceStabilizingTransformation(dds, blind=FALSE)
res$miRNA <- row.names(res)
col_order <- c("miRNA","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
res2 <- res[,col_order]
write.table(res2, file=args[5], sep="\t", row.names = FALSE)
#write.table(res, "deseq_output.txt", sep="\t")

pdf(file=args[4], width=8,height=8)
#pdf("pca_plot.pdf", width=8,height=8)
plotPCA(vsdata, intgroup="group")
garbage <- dev.off()
save.image(file=args[6]) 
