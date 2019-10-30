

# 1. Load all libraries ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximportData")
BiocManager::install("tximport")

library("tidyverse")
library("tximport")
library("readr")
library("tximportData")

BiocManager::install("DESeq2")
library("DESeq2")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

BiocManager::install("DRIMSeq")
library("DRIMSeq")
# 2. Read in files --------------------------------------------------------

setwd("~/Desktop/ciona_robusta_2019_sept/ciona_zf/")
dir <- "data/"
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE, fill = TRUE)
samples_removeRep2 <- read.table(file.path(dir,"samples_minus_rep2.txt"), header=TRUE, fill = TRUE)
samples_switchRep2 <- read.table(file.path(dir,"samples.txt"), header=TRUE, fill = TRUE)
sample1 <- read.table(file.path(dir,"samples1.txt"), header=TRUE, fill = TRUE)
sample2 <- read.table(file.path(dir,"samples2.txt"), header=TRUE, fill = TRUE)
sample3 <- read.table(file.path(dir,"samples3.txt"), header=TRUE, fill = TRUE)

samples$condition <- factor(rep(c("WT","Zf_cpr"),each=12))
samples_removeRep2$condition <- factor(rep(c("WT","Zf_cpr"),each=8))
samples_switchRep2$condition <- factor(rep(c("WT","Zf_cpr"),each=4))
sample1$condition <- factor(rep(c("WT","Zf_cpr"),each=4))
sample2$condition <- factor(rep(c("WT","Zf_cpr"),each=4))
sample3$condition <- factor(rep(c("WT","Zf_cpr"),each=4))

rownames(samples) <- samples$assay
rownames(samples_removeRep2) <- samples_removeRep2$assay
rownames(samples_switchRep2) <- samples_switchRep2$assay
rownames(sample1) <- sample1$assay
rownames(sample2) <- sample2$assay
rownames(sample3) <- sample3$assay

sample3[,c("assay","sample","rep","condition")]
samples_removeRep2[,c("assay","sample","rep","condition")]
samples_switchRep2[,c("assay","sample","rep","condition")]

files <- file.path(dir,"quants", paste0(samples$assay,"_quant"), "quant.sf")
files_removeRep2 <- file.path(dir,"quants", paste0(samples_removeRep2$assay,"_quant"), "quant.sf")
files_switchRep2 <- file.path(dir,"quants", paste0(samples_switchRep2$assay,"_quant"), "quant.sf")
files1 <- file.path(dir,"quants", paste0(sample1$assay,"_quant"), "quant.sf")
files2 <- file.path(dir,"quants", paste0(sample2$assay,"_quant"), "quant.sf")
files3 <- file.path(dir,"quants", paste0(sample3$assay,"_quant"), "quant.sf")

names(files) <- samples$assay
names(files_removeRep2) <- samples_removeRep2$assay
names(files_switchRep2) <- samples_switchRep2$assay
names(files1) <- sample1$assay
names(files2) <- sample2$assay
names(files3) <- sample3$assay
all(file.exists(files))
#tx2geneKH <- data.frame(TXNAME = rownames(dds), GENEID = gsub("\\.v.+", "", rownames(dds)))
#tx2geneNCBI <- read_csv("NCBITx2gene.txt", col_names = c("TXNAME","GENEID"))
#tx2gene <- read_csv("data/ciona_rob.tx2gene2.csv", col_names = c("GENEID","TXNAME")) %>% dplyr::select(TXNAME, GENEID)
tx2gene <- read_csv("data/ciona_rob.tx2gene.csv")
txi.g <- tximport(files, type="salmon", tx2gene=tx2gene)
txi_removeRep2.g <- tximport(files_removeRep2, type="salmon", tx2gene=tx2gene)
txi_switchRep2.g <- tximport(files_switchRep2, type="salmon", tx2gene=tx2gene)
txi1.g <- tximport(files1, type="salmon", tx2gene=tx2gene)
txi2.g <- tximport(files2, type="salmon", tx2gene=tx2gene)
txi3.g <- tximport(files3, type="salmon", tx2gene=tx2gene)
txi.t <- tximport(files, type="salmon", txOut = TRUE, countsFromAbundance="scaledTPM")
txi_removeRep2.t <- tximport(files_removeRep2, type="salmon", txOut = TRUE, countsFromAbundance="scaledTPM")
txi_switchRep2.t <- tximport(files_switchRep2, type="salmon", txOut = TRUE, countsFromAbundance="scaledTPM")
names(txi)
#write_csv(tx2gene, "data/ciona_rob.tx2gene.csv")

# 3. Setup DE dataframes --------------------------------------------------

dds.g <- DESeqDataSetFromTximport(txi.g, samples, ~ rep + condition)
dds_removeRep2.g <- DESeqDataSetFromTximport(txi_removeRep2.g, samples_removeRep2, ~ rep + condition)
dds_switchRep2.g <- DESeqDataSetFromTximport(txi_switchRep2.g, samples_switchRep2, ~ rep + condition)
dds1.g <- DESeqDataSetFromTximport(txi1.g, sample1, ~ lane + condition)
dds2.g <- DESeqDataSetFromTximport(txi2.g, sample2, ~ lane + condition)
dds3.g <- DESeqDataSetFromTximport(txi3.g, sample3, ~ lane + condition)

dds <- collapseReplicates(dds.g, dds.g$sample)
dds_removeRep2.g <- collapseReplicates(dds_removeRep2.g, dds_removeRep2.g$sample)
dds_switchRep2.g <- collapseReplicates(dds_switchRep2.g, dds_switchRep2.g$sample)
dds1.g <- collapseReplicates(dds1.g, dds1.g$sample)
dds2.g <- collapseReplicates(dds2.g, dds2.g$sample)
dds3.g <- collapseReplicates(dds3.g, dds3.g$sample)

boxplot(log10(counts(dds)+1))
dds <- estimateSizeFactors(dds)
boxplot(log10(counts(dds,normalized=TRUE)+1))

vsd.g <- vst(dds.g)
pcaData <- plotPCA(vsd.g, intgroup=c("condition", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
svg("all_gene-level_pca.svg")
#plotPCA(vsd.g, c("condition","rep"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
mat.g <- assay(vsd.g)
mat.g <- limma::removeBatchEffect(mat.g, vsd.g$rep)
assay(vsd.g) <- mat.g
svg("all_gene-level_pca_batch.svg")
plotPCA(vsd.g, c("condition","rep"))
dev.off()

vsd_removeRep2.g <- vst(dds_removeRep2.g)
pcaData <- plotPCA(vsd_removeRep2.g, intgroup=c("condition", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
svg("removeRep2_gene-level_pca.svg")
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
plotPCA(vsd_removeRep2.g, c("condition","rep"))
mat_removeRep2.g <- assay(vsd_removeRep2.g)
mat_removeRep2.g <- limma::removeBatchEffect(mat_removeRep2.g, vsd_removeRep2.g$rep)
assay(vsd_removeRep2.g) <- mat_removeRep2.g
pcaData <- plotPCA(vsd_removeRep2.g, intgroup=c("condition", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
svg("removeRep2_gene-level_pca_batch.svg")
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
plotPCA(vsd_removeRep2.g)

vsd_switchRep2.g <- vst(dds_switchRep2.g)
pcaData <- plotPCA(vsd_switchRep2.g, intgroup=c("condition", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
svg("switchRep2_gene-level_pca.svg")
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
plotPCA(vsd_switchRep2.g, c("condition","rep"))
mat_switchRep2.g <- assay(vsd_switchRep2.g)
mat_switchRep2.g <- limma::removeBatchEffect(mat_switchRep2.g, vsd_switchRep2.g$rep)
assay(vsd_switchRep2.g) <- mat_switchRep2.g
pcaData <- plotPCA(vsd_switchRep2.g, intgroup=c("condition", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
svg("switchRep2_gene-level_pca_batch.svg")
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
plotPCA(vsd_switchRep2.g)

dds.g$condition <- relevel(dds.g$condition, ref = "WT")
dds.g <- DESeq(dds.g)
res.g <- DESeq2::results(dds.g)
#res.g
summary(res.g)
head(res.g[order(res.g$pvalue),])

dds1.g$condition <- relevel(dds1.g$condition, ref = "WT")
dds1.g <- DESeq(dds1.g)
res1.g <- DESeq2::results(dds1.g)
#res.g
summary(res1.g)
head(res1.g[order(res1.g$pvalue),])

dds2.g$condition <- relevel(dds2.g$condition, ref = "WT")
dds2.g <- DESeq(dds2.g)
res2.g <- DESeq2::results(dds2.g)
#res.g
summary(res2.g)
head(res2.g[order(res2.g$pvalue),])

dds3.g$condition <- relevel(dds3.g$condition, ref = "WT")
dds3.g <- DESeq(dds3.g)
res3.g <- DESeq2::results(dds3.g)
#res.g
summary(res3.g)
head(res3.g[order(res3.g$pvalue),])

res1.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res1.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange > 0)
  
rep1.up <- res1.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res1.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange > 0)
rep2.up <- res2.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res2.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange > 0)
rep3.up <- res3.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res3.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange > 0)

rep1.down <- res1.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res1.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange < 0)
rep2.down <- res2.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res2.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange < 0)
rep3.down <- res3.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res3.g)) %>% 
  filter(pvalue < 0.05 & log2FoldChange < 0)

listInput <- list("Rep 1 up" = rep1.up$Kh.id, "Rep 1 down" = rep1.down$Kh.id,
                  "Rep 2 up" = rep2.up$Kh.id, "Rep 2 down" = rep2.down$Kh.id,
                  "Rep 3 up" = rep3.up$Kh.id, "Rep 3 down" = rep3.down$Kh.id)

dds_removeRep2.g$condition <- relevel(dds_removeRep2.g$condition, ref = "WT")
dds_removeRep2.g <- DESeq(dds_removeRep2.g)
res_removeRep2.g <- DESeq2::results(dds_removeRep2.g)
#res_removeRep2.g
summary(res_removeRep2.g)
head(res_removeRep2.g[order(res_removeRep2.g$pvalue),])

dds_switchRep2.g$condition <- relevel(dds_switchRep2.g$condition, ref = "WT")
dds_switchRep2.g <- DESeq(dds_switchRep2.g)
res_switchRep2.g <- DESeq2::results(dds_switchRep2.g)
#res
summary(res_switchRep2.g)
head(res_switchRep2.g[order(res_switchRep2.g$pvalue),])
plotMA(res_switchRep2.g)

resultsNames(dds)
resApeT <- lfcShrink(dds, coef="condition_Zf_cpr_vs_WT", type="apeglm")
plotMA(resApeT)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
summary(resApeT)


# Volcano plots -----------------------------------------------------------

svg("gene-level_vol.svg")
EnhancedVolcano(res1.g,
                lab = rownames(res.g),
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()

svg("gene-level_vol_removeR2.svg")
EnhancedVolcano(res_removeRep2.g,
                lab = rownames(res_removeRep2.g),
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()

svg("gene-level_vol_switchR2.svg")
EnhancedVolcano(res_switchRep2.g,
                lab = rownames(res_switchRep2.g),
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()


# Combine with in-situ data -----------------------------------------------


library(readxl)
insitu_data <- read_excel("~/Desktop/ciona_robusta_2019_sept/ciona_zf/data/KHID-UniqueName-URLs-InSitu-COMPLETE.xlsx") #readxl library

resSorted <- res_switchRep2.g[order(res_switchRep2.g$pvalue),]
resSorted %>%
  as_tibble() %>% 
  mutate(KHID = rownames(resSorted)) %>% 
  dplyr::select(KHID, baseMean, log2FoldChange, pvalue, padj) %>% 
  left_join(insitu_data) -> resISH
resISH %>% View()

write_csv(resISH, "Crob_zfp36_crispr2.csv")

res.g_0.05 <- res.g %>% as_tibble() %>% mutate(KHid = rownames(res.g)) %>% filter(pvalue <= 0.05)
res_removeRep2.g_0.05 <- res_removeRep2.g %>% as_tibble() %>% mutate(KHid = rownames(res_removeRep2.g)) %>% filter(pvalue <= 0.05)
res_switchRep2.g_0.05 <- res_switchRep2.g %>% as_tibble() %>% mutate(KHid = rownames(res_switchRep2.g)) %>% filter(pvalue <= 0.05)

sum(res.g_0.05$KHid %in% res_removeRep2.g_0.05$KHid & res.g_0.05$KHid %in% res_switchRep2.g_0.05$KHid)
sum(res.g_0.05$KHid %in% res_removeRep2.g_0.05$KHid & !(res.g_0.05$KHid %in% res_switchRep2.g_0.05$KHid))
sum(!(res.g_0.05$KHid %in% res_removeRep2.g_0.05$KHid) & (res.g_0.05$KHid %in% res_switchRep2.g_0.05$KHid))
sum(!(res.g_0.05$KHid %in% res_removeRep2.g_0.05$KHid) & !(res.g_0.05$KHid %in% res_switchRep2.g_0.05$KHid))
sum(!(res_switchRep2.g_0.05$KHid %in% res_removeRep2.g_0.05$KHid) & !(res_switchRep2.g_0.05$KHid %in% res.g_0.05$KHid))
sum(!(res_removeRep2.g_0.05$KHid %in% res_switchRep2.g_0.05$KHid) & !(res_removeRep2.g_0.05$KHid %in% res.g_0.05$KHid))
sum((res_switchRep2.g_0.05$KHid %in% res_removeRep2.g_0.05$KHid) & !(res_switchRep2.g_0.05$KHid %in% res.g_0.05$KHid))

res.g_0.05a <- res.g %>% as_tibble() %>% mutate(KHid = rownames(res.g)) %>% filter(padj <= 0.05)
res_removeRep2.g_0.05a <- res_removeRep2.g %>% as_tibble() %>% mutate(KHid = rownames(res_removeRep2.g)) %>% filter(padj <= 0.05)
res_switchRep2.g_0.05a <- res_switchRep2.g %>% as_tibble() %>% mutate(KHid = rownames(res_switchRep2.g)) %>% filter(padj <= 0.05, log2FoldChange < -1)

sum(res.g_0.05a$KHid %in% res_removeRep2.g_0.05a$KHid & res.g_0.05a$KHid %in% res_switchRep2.g_0.05a$KHid)
sum(res.g_0.05a$KHid %in% res_removeRep2.g_0.05a$KHid & !(res.g_0.05a$KHid %in% res_switchRep2.g_0.05a$KHid))
sum(!(res.g_0.05a$KHid %in% res_removeRep2.g_0.05a$KHid) & (res.g_0.05a$KHid %in% res_switchRep2.g_0.05a$KHid))
sum(!(res.g_0.05a$KHid %in% res_removeRep2.g_0.05a$KHid) & !(res.g_0.05a$KHid %in% res_switchRep2.g_0.05a$KHid))
sum(!(res_switchRep2.g_0.05a$KHid %in% res_removeRep2.g_0.05a$KHid) & !(res_switchRep2.g_0.05a$KHid %in% res.g_0.05a$KHid))
sum(!(res_removeRep2.g_0.05a$KHid %in% res_switchRep2.g_0.05a$KHid) & !(res_removeRep2.g_0.05a$KHid %in% res.g_0.05a$KHid))
sum((res_switchRep2.g_0.05a$KHid %in% res_removeRep2.g_0.05a$KHid) & !(res_switchRep2.g_0.05a$KHid %in% res.g_0.05a$KHid))
