#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
print("Program path: /programs/cMYC-induced.miR.genes.comaprison.r")
if(length(args) != 1) {stop("The path to the EDGER Robject is missing.")}
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
#Robject.path <- '/workspace/lukasz/Michalina_Z/RWorkspace_miRNA/BL_BLL.RData'
condition <- sub(sub(Robject.path, pattern = "^.*RWorkspace_miRNA/", replacement = ""), pattern = "\\.RData$", replacement = "")
workdir <- paste("/workspace/lukasz/Michalina_Z/New_analyses/miR", condition, sep = "/")
dir.create(workdir, recursive = T)
print(workdir)
setwd(workdir)
load(Robject.path)
save.image("cMYC.induced.miR.genes.RData")
con <- file('/workspace/lukasz/Michalina_Z/New_analyses/BL_SIGNATURE.miR.txt')
gene.list1 <- readLines(con)
close(con)
gene.list1 <- gene.list1[!duplicated(gene.list1)]

gene.list1.conv <- NULL
for(i in gene.list1) {if(length(alias2Symbol(i, species = "Hs")) == 0) {gene.list1.conv <- append(gene.list1.conv, values = i)} else {gene.list1.conv <- append(gene.list1.conv, values = alias2Symbol(i, species = "Hs"))}}

colnames(final_columns) <- gsub(gsub(gsub(colnames(final_columns), pattern = "^[^0-9]+", replacement = ""), pattern = "[^0-9]+$", replacement = ""), pattern = "i?_[0-9]+", replacement = "")
colnames(final_columns)[colnames(final_columns) == "126"] <- "1263"
metadata$SampleID[metadata$SampleID == "126"] <- "1263"

del.cols <- c("186", "140", "227")
del.cols <- del.cols[del.cols %in% colnames(final_columns)]
if(length(del.cols) > 0) { final_columns <- final_columns %>% dplyr::select(!del.cols)}

samples.list <- as.data.frame(colnames(final_columns))
colnames(samples.list) <- c("SampleID")
metadata[["SampleID"]] <- gsub(gsub(gsub(metadata[["SampleID"]], pattern = "^[^0-9]+", replacement = ""), pattern = "[^0-9]+$", replacement = ""), pattern = "i?_[0-9]+", replacement = "")

sampleTable2 <- merge(x = samples.list, y = metadata, by.x = "SampleID", by.y = "SampleID")
rownames(sampleTable2) <- sampleTable2$SampleID
design <- model.matrix(~sampleTable2[["condition"]])

if (length(sampleTable2[["SampleID"]])<3) {
  save.image(file=paste0(condition,".SE",".RData"))
  stop('The analysis cannot be performed for less than three samples.')}

library(edgeR)
y <- DGEList(counts=final_columns)
# keep <- rowSums(cpm(y)>1) >= 2
# y <- y[keep, , keep.lib.sizes=FALSE]
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
rld <- cpm(y, log=TRUE)
edgeR.miR.counts.normalized <- cpm(y, log=FALSE)

rld.sel <- rld[rownames(rld) %in% gene.list1.conv,]
rld.sel <- rld.sel[!duplicated(rownames(rld.sel)),]

anno <- sampleTable2 %>% dplyr::select("condition")
pdf(paste(condition, "cMYC-induced miR genes - BL_signature comparison heatmap.pdf", sep = "-"))
heatmap0 <- pheatmap(rld.sel, annotation_col = anno)
dev.off()

rld.sel.medians <- rowMedians(rld.sel)
names(rld.sel.medians) <- rownames(rld.sel)
rld.sel.cat <- rld.sel
for(i in seq(1,length(rld.sel.medians))) {j <- rld.sel.medians[i]
rld.sel.cat[names(j),][rld.sel[names(j),] > j] <- 1
rld.sel.cat[names(j),][rld.sel[names(j),] <= j] <- 0
}
rld.sel.cat.t <- t(rld.sel.cat)
rld.sel.cat.t <- rld.sel.cat.t[order(rownames(rld.sel.cat.t)),]
anno <- anno %>% dplyr::arrange(rownames(anno))

if(all(rownames(rld.sel.cat.t) == rownames(anno)) == TRUE) {
  rld.sel.cat.t.by <- by(data = rld.sel.cat.t, INDICES = anno$condition, FUN = colMeans)} else 
  { save.image("cMYC.induced.miR.genes.RData")
    stop("Annotations do not match sample names.")}
rld.sel.cat.t.comb <- foreach(name = names(rld.sel.cat.t.by), .combine = "cbind") %do% {rld.sel.cat.t.by[[name]]}
colnames(rld.sel.cat.t.comb) <- names(rld.sel.cat.t.by)
heatmap1 <- pheatmap(rld.sel.cat.t.comb, main = "cMYC-induced miR genes (median expression-categorized, condition-grouped)", filename = paste(condition, "cMYC-induced miR genes (median expression-categorized, condition-grouped).BL_signature.heatmap.pdf", sep = "."))
rld.sel.cat.t.comb.t <- t(rld.sel.cat.t.comb)
rld.sel.cat.t.comb.t.pca.data <- prcomp(rld.sel.cat.t.comb.t)
rld.sel.cat.t.comb.t.df <- as.data.frame(rld.sel.cat.t.comb.t)
rld.sel.cat.t.comb.t.df$Condition <- rownames(rld.sel.cat.t.comb.t.df)

pdf(paste(condition, "cMYC-induced miR genes (median expression-categorized, condition-grouped).BL_signature.PCA.pdf", sep = "."))
autoplot(rld.sel.cat.t.comb.t.pca.data, data = rld.sel.cat.t.comb.t.df, size = 3, colour = "Condition") + 
  scale_color_manual(values = if (nrow(rld.sel.cat.t.comb.t.pca.data$x)>25) {rainbow(nrow(rld.sel.cat.t.comb.t.pca.data$x))} else {as.vector(cols25())}) + coord_fixed() + 
  labs(title = "PCA plot for cMYC-induced miR genes") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

rld.sel.cat.t.comb[rld.sel.cat.t.comb >= 0.5] <- 1
rld.sel.cat.t.comb[rld.sel.cat.t.comb < 0.5] <- 0
rld.sel.cat.t.comb.df <- as.data.frame(rld.sel.cat.t.comb) 
heatmap2 <- pheatmap(rld.sel.cat.t.comb, main = "cMYC-induced miR genes (median expression-categorized, binarized, condition-grouped)", filename = paste(condition, "cMYC-induced miR genes (median expression-categorized, binarized, condition-grouped).BL_signature.heatmap.pdf", sep = "."))

if("BL" %in% colnames(rld.sel.cat.t.comb.df)) {
rld.sel.cat.t.comb.df.bool <- rld.sel.cat.t.comb.df == rld.sel.cat.t.comb.df[["BL"]]
rld.sel.cat.t.comb.df.bool.t <- t(rld.sel.cat.t.comb.df.bool)

f_fraction <- function(x) {
  return(paste(rowSums(x), "matching out of", length(x), "genes", "(", rowSums(x)/length(x), ")"))}

fin.summary <- by(rld.sel.cat.t.comb.df.bool.t, INDICES = rownames(rld.sel.cat.t.comb.df.bool.t), FUN = f_fraction)
sink(file =  paste(condition, "cMYC-induced miR genes-summary.BL_signature.txt", sep = "."))
print(fin.summary)
sink() }

j <- NULL
for(i in seq(1, nrow(rld.sel))){if(all(rld.sel[i,] == 0)){j <- append(j, values = i)}}
if(!is.null(j)) {rld.sel.no.zeros <- rld.sel[-j,]} else {rld.sel.no.zeros <- rld.sel}
rld.sel.no.zeros.txp <- t(rld.sel.no.zeros)
rld.sel.no.zeros.txp.merged <- merge(x = rld.sel.no.zeros.txp, y = sampleTable2, by.x = 0, by.y = 0)
delcols <- c(1, ncol(rld.sel.no.zeros.txp.merged), ncol(rld.sel.no.zeros.txp.merged)-1)

pca.data <- prcomp(rld.sel.no.zeros.txp.merged[,-delcols])

pdf(paste(condition, "cMYC-induced miR genes - BL_signature PCA plot.pdf", sep = "-"))
autoplot(pca.data, data = rld.sel.no.zeros.txp.merged, shape = "condition", colour = "SampleID", size = 3) + 
  scale_color_manual(values = if (nrow(pca.data$x)>25) {rainbow(nrow(pca.data$x))} else {as.vector(cols25())}) + coord_fixed() + 
  labs(title = "PCA plot for cMYC-induced miR genes") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
save.image("cMYC.induced.miR.genes.RData")

cat("The analysis is complete.\n")