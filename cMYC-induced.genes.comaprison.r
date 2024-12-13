#! /usr/bin/env Rscript

cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 1) {stop("The path to the DESEQ2 Robject is missing.")}
library(limma)
library(foreach)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(dplyr)
library(pals)
library(ggfortify)
Robject.path <- args[1]
path.elements <- unlist(strsplit(Robject.path, split = "/"))
workdir <- paste(path.elements[seq(1,length(path.elements)-1)], collapse = "/")
print(workdir)
setwd(workdir)
condition <- sub(sub(Robject.path, pattern = "^.*ALL_BAMS/", replacement = ""), pattern = "\\/DESEQ2\\/.*$", replacement = "")
Robject.path <- path.elements[length(path.elements)]

load(Robject.path)
signature.path <- '/workspace/lukasz/NGS-all-in-one/TXT_FILES/MYC_Seitz.2011.txt'
signature.name <- sub(basename(signature.path), pattern = "\\.txt", replacement = "")
con <- file(signature.path)
gene.list1 <- readLines(con)
close(con)

gene.list1.conv <- NULL
for(i in gene.list1) {if(length(alias2Symbol(i, species = "Hs")) == 0) {gene.list1.conv <- append(gene.list1.conv, values = i)} else {gene.list1.conv <- append(gene.list1.conv, values = alias2Symbol(i, species = "Hs"))}}

rld.mx <- assay(rld)
rownames(rld.mx) <- sub(rownames(rld.mx), pattern = "\\.[0-9]+$", replacement = "")

rownames(rld.mx) <- as.vector(mapIds(org.Hs.eg.db,
                 keys=rownames(rld.mx),
                 column="SYMBOL",
                 keytype="ENSEMBL",
                 multiVals="first"))

rld.mx <- rld.mx[!is.na(rownames(rld.mx)),]
rld.mx.sel <- rld.mx[rownames(rld.mx) %in% gene.list1.conv,]
rld.mx.sel <- rld.mx.sel[!duplicated(rownames(rld.mx.sel)),]

anno <- sampleTable %>% dplyr::select("condition")
pdf(paste(condition, "cMYC-induced genes", signature.name, "comparison heatmap.pdf", sep = "-"), height = 15, width = 8)
heatmap0 <- pheatmap(rld.mx.sel, annotation_col = anno, fontsize_row = 5)
dev.off()

rld.mx.sel.medians <- rowMedians(rld.mx.sel)
names(rld.mx.sel.medians) <- rownames(rld.mx.sel)
rld.mx.sel.cat <- rld.mx.sel
for(i in seq(1,length(rld.mx.sel.medians))) {j <- rld.mx.sel.medians[i]
rld.mx.sel.cat[names(j),][rld.mx.sel[names(j),] > j] <- 1
rld.mx.sel.cat[names(j),][rld.mx.sel[names(j),] <= j] <- 0
}
rld.mx.sel.cat.t <- t(rld.mx.sel.cat)
if(all(rownames(rld.mx.sel.cat.t) == rownames(anno)) == TRUE) {
  rld.mx.sel.cat.t.by <- by(data = rld.mx.sel.cat.t, INDICES = anno$condition, FUN = colMeans)} else 
  {stop("Annotations do not match sample names.")}
rld.mx.sel.cat.t.comb <- foreach(name = names(rld.mx.sel.cat.t.by), .combine = "cbind") %do% {rld.mx.sel.cat.t.by[[name]]}
colnames(rld.mx.sel.cat.t.comb) <- names(rld.mx.sel.cat.t.by)
heatmap1 <- pheatmap(rld.mx.sel.cat.t.comb, main = "cMYC-induced genes (median expression-categorized, condition-grouped)", filename = paste(condition, "cMYC-induced genes (median expression-categorized, condition-grouped)", signature.name, "heatmap.pdf", sep = "."), height = 15, width = 8, fontsize_row = 5)
rld.mx.sel.cat.t.comb.t <- t(rld.mx.sel.cat.t.comb)
rld.mx.sel.cat.t.comb.t.pca.data <- prcomp(rld.mx.sel.cat.t.comb.t)
rld.mx.sel.cat.t.comb.t.df <- as.data.frame(rld.mx.sel.cat.t.comb.t)
rld.mx.sel.cat.t.comb.t.df$Condition <- rownames(rld.mx.sel.cat.t.comb.t.df)

pdf(paste(condition, "cMYC-induced genes (median expression-categorized, condition-grouped)", signature.name, "PCA.pdf", sep = "."))
autoplot(rld.mx.sel.cat.t.comb.t.pca.data, data = rld.mx.sel.cat.t.comb.t.df, size = 3, colour = "Condition") + 
  scale_color_manual(values = if (nrow(rld.mx.sel.cat.t.comb.t.pca.data$x)>25) {rainbow(nrow(rld.mx.sel.cat.t.comb.t.pca.data$x))} else {as.vector(cols25())}) + coord_fixed() + 
  labs(title = "PCA plot for cMYC-induced genes") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

rld.mx.sel.cat.t.comb[rld.mx.sel.cat.t.comb >= 0.5] <- 1
rld.mx.sel.cat.t.comb[rld.mx.sel.cat.t.comb < 0.5] <- 0
rld.mx.sel.cat.t.comb.df <- as.data.frame(rld.mx.sel.cat.t.comb) 
heatmap2 <- pheatmap(rld.mx.sel.cat.t.comb, main = "cMYC-induced genes (median expression-categorized, binarized, condition-grouped)", filename = paste(condition, "cMYC-induced genes (median expression-categorized, binarized, condition-grouped)", signature.name, "heatmap.pdf", sep = "."), height = 15, width = 8, fontsize_row = 5)
rld.mx.sel.cat.t.comb.df.bool <- rld.mx.sel.cat.t.comb.df == rld.mx.sel.cat.t.comb.df$BL
rld.mx.sel.cat.t.comb.df.bool.t <- t(rld.mx.sel.cat.t.comb.df.bool)

f_fraction <- function(x) {
  return(paste(rowSums(x), "matching out of", length(x), "genes", "(", rowSums(x)/length(x), ")"))}

fin.summary <- by(rld.mx.sel.cat.t.comb.df.bool.t, INDICES = rownames(rld.mx.sel.cat.t.comb.df.bool.t), FUN = f_fraction)
sink(file =  paste(condition, "cMYC-induced genes-summary", signature.name, "txt", sep = "."))
print(fin.summary)
sink()

j <- NULL
for(i in seq(1, nrow(rld.mx.sel))){if(all(rld.mx.sel[i,] == 0)){j <- append(j, values = i)}}
rld.mx.sel.no.zeros <- rld.mx.sel[-j,]
rld.mx.sel.no.zeros.txp <- t(rld.mx.sel.no.zeros)
rld.mx.sel.no.zeros.txp.merged <- merge(x = rld.mx.sel.no.zeros.txp, y = sampleTable, by.x = 0, by.y = 0)
delcols <- c(1, ncol(rld.mx.sel.no.zeros.txp.merged), ncol(rld.mx.sel.no.zeros.txp.merged)-1)

pca.data <- prcomp(rld.mx.sel.no.zeros.txp.merged[,-delcols])

pdf(paste(condition, "cMYC-induced genes", signature.name, "PCA plot.pdf", sep = "-"))
autoplot(pca.data, data = rld.mx.sel.no.zeros.txp.merged, shape = "condition", colour = "SampleID", size = 3) + 
  scale_color_manual(values = if (nrow(pca.data$x)>25) {rainbow(nrow(pca.data$x))} else {as.vector(cols25())}) + coord_fixed() + 
  labs(title = "PCA plot for cMYC-induced genes") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
save.image("cMYC.induced.genes.RData")

cat("The analysis is complete.\n")