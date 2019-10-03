

# 1. Load all libraries ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximportData")
BiocManager::install("tximport")

library("tximport")
library("readr")
library("tximportData")

BiocManager::install("DESeq2")
library("DESeq2")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

# 2. Read in files --------------------------------------------------------

setwd("~/Desktop/ciona_robusta_2019_sept/ciona_zf/")
dir <- "data/"
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE, fill = TRUE)
samples$condition <- factor(rep(c("WT","Zf_cpr"),each=12))
rownames(samples) <- samples$assay #run
samples[,c("assay","sample","rep","condition")]

files <- file.path(dir,"quants", paste0(samples$assay,"_quant"), "quant.sf")
names(files) <- samples$assay
all(file.exists(files))
tx2gene <- data.frame(TXNAME = rownames(dds), GENEID = gsub("\\.v.+", "", rownames(dds)))
tx2gene <- read_csv("data/ciona_rob.tx2gene.cscv")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

write_csv(tx2gene, "ciona_rob.tx2gene.cscv")

dds <- DESeqDataSetFromTximport(txi, samples, ~ rep + condition)
dds <- collapseReplicates(dds, dds$sample)

keep <- rowSums(counts(dds1) >= 5) >= 3
table(keep)
dds <- dds[keep,]

boxplot(log10(counts(dds1)+1))
dds <- estimateSizeFactors(dds1)
boxplot(log10(counts(dds,normalized=TRUE)+1))

vsd <- vst(dds)
plotPCA(vsd, c("condition","rep"))
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$rep)
assay(vsd) <- mat
plotPCA(vsd)

dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
head(res[order(res$pvalue),])

plotMA(res)

resultsNames(dds)
resApeT <- lfcShrink(dds, coef="condition_Zf_cpr_vs_WT", type="apeglm")
plotMA(resApeT)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
summary(resApeT)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(res,
                lab = rownames(resApeT),
                x = 'log2FoldChange',
                y = 'pvalue')

library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","rep")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

