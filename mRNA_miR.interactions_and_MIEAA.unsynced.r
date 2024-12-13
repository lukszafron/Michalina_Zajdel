#! /usr/bin/env Rscript
# for i in $(find -iregex "^\.\/.+\.RData$"); do NAME=$(basename $i);echo -e ""$NAME" is being analyzed...\n" ; R CMD BATCH --no-restore --no-save "--args $NAME" DESeq2.R $NAME.log; if [[ $? -eq 0 ]]; then echo "$NAME is complete"; else echo "$NAME was terminated with an error."; fi; done
print("Program path: /programs/mRNA_miR.interactions_and_MIEAA.r")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {stop("Two arguments are obligatory. The first should point to the RData file with the final DESeq2 expression results. The second should show the RData file with the final miRNA results obtained by M. Kulecka.")}

mrnafile <- args[1]
path.elements <- unlist(strsplit(mrnafile, split = "/"))
workdir <- paste(path.elements[seq(1,length(path.elements)-1)], collapse = "/")
dir.create(paste(workdir, "miRNA", sep = "/"))
setwd(paste(workdir, "miRNA", sep = "/"))
condition <- sub(sub(mrnafile, pattern = "^.*ALL_BAMS/", replacement = ""), pattern = "/DESEQ2/UNFILTERED.*", replacement = "")
mrnafile <- path.elements[length(path.elements)]
mirfile <- args[2]

#Filters <- "UNFILTERED"
#Filters <- "FILTERED"
#workspace <- '/workspace/lukasz/Michalina Z/New_analyses/BL_BLL/BAMS/DESEQ2 - without 3 samples, unfiltered/DESEQ2.SE.transformed.filtered.RData'
#mirdir <- '/workspace/lukasz/Michalina Z/RWorkspace_miRNA'
load(paste("..",mrnafile, sep = "/")) #mRNA RData
load(mirfile) # miR RData

colnames(final_columns) <- gsub(gsub(gsub(colnames(final_columns), pattern = "^[^0-9]+", replacement = ""), pattern = "[^0-9]+$", replacement = ""), pattern = "i?_[0-9]+", replacement = "")
colnames(final_columns)[colnames(final_columns) == "126"] <- "1263"
metadata$SampleID[metadata$SampleID == "126"] <- "1263"
library(dplyr)
del.cols <- c("186", "140", "227")
del.cols <- del.cols[del.cols %in% colnames(final_columns)]
if(length(del.cols) > 0) { final_columns <- final_columns %>% dplyr::select(!del.cols)}

samples.list <- as.data.frame(colnames(final_columns))
colnames(samples.list) <- c(sample_col)
metadata[[sample_col]] <- gsub(gsub(gsub(metadata[[sample_col]], pattern = "^[^0-9]+", replacement = ""), pattern = "[^0-9]+$", replacement = ""), pattern = "i?_[0-9]+", replacement = "")

sampleTable2 <- merge(x = samples.list, y = metadata, by.x = sample_col, by.y = "SampleID")

design <- model.matrix(~sampleTable2[[ind.factor]])

if (length(sampleTable2[[sample_col]])<3) {
  save.image(file=paste0(condition,".SE.transformed",".RData"))
  stop('The analysis cannot be performed for less than three samples.')}

library(edgeR)
y <- DGEList(counts=final_columns)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y, design)
#Testowanie, czy są jakieś różnice między grupami, bez definiowania, o jakie grupy chodzi
lrt <- glmLRT(fit)
data.miR<-lrt$table[with(lrt$table,order(PValue)),]
data.miR$FDR<-p.adjust(data.miR[["PValue"]],method="fdr")
data.miR$FC <- 2**data.miR$logFC
data.miR <- data.miR[with(data.miR,order(FDR)),]
data.miR$miR.names <- rownames(data.miR)
data.miR <- data.miR[,c(7, 2:6)]
data.miR.05 <- subset(data.miR, subset = FDR < 0.05)

write.table(data.miR, file = paste("edgeR_analysis_results",".(",condition,").","miR",".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
library("ReportingTools")
htmlRep <- HTMLReport(shortName="edgeR analysis report", title=paste("edgeR_analysis_results",".(",condition,").","miR", sep = ""),
reportDirectory=paste("edgeR_analysis_results",".(",condition,").","miR", sep = ""))

publish(data.miR, htmlRep)
url <- finish(htmlRep)
browseURL(url, browser="firefox")

rld <- cpm(y, log=TRUE)
edgeR.miR.counts.normalized <- cpm(y, log=FALSE)

  pdf(paste("Read counts normalization-miR",condition,".pdf", sep = "."))
  par(mfrow=c(2,1))
  boxplot(y$counts+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
  title(main=paste("edgeR - raw read counts - miR", condition, sep = "."), xlab=sample_col, ylab="Read counts in a log10 scale")
  boxplot(edgeR.miR.counts.normalized+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
  title(main=paste("edgeR - normalized read counts - miR", condition, sep = "." ), xlab=sample_col, ylab="Read counts in a log10 scale")
  dev.off()

library(pals)
library(ggplot2)
plot.mds <- plotMDS(y)
dev.off()
unlink('Rplots.pdf')
plot.mds.df <- as.data.frame(plot.mds$cmdscale.out)
  pdf(paste("MDS-plot",".(",condition,").","edgeR",".pdf", sep = ""))
  if(nrow(plot.mds.df)<=25) {
    ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = sampleTable2[[sample_col]], shape = sampleTable2[["condition"]])) + 
      geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
      xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
      labs(title = "Multidimensional scaling (MDS) plot", shape = condition, color = sample_col) + 
      theme(plot.title = element_text(hjust = 0.5))
    print(ggplot.1)
  } else 	{
    ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = sampleTable2[["condition"]])) + 
      geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
      xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
      labs(title = "Multidimensional scaling (MDS) plot", color = condition) + 
      theme(plot.title = element_text(hjust = 0.5))
    print(ggplot.1)
  }
  dev.off()

  library("genefilter")
  #library(tidyverse)
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  
  topVarGenes <- head(order(rowVars(rld), decreasing = TRUE), 20)
  mat  <- rld[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
 
  library("pheatmap")
  library("RColorBrewer")

  anno <- as.data.frame(sampleTable2[, ind.factor],row.names = as.vector(sampleTable2[[sample_col]]))
  colnames(anno) <- ind.factor
  pdf(paste("TopVarGenes.rld",".(",condition,").","miR",".pdf", sep = ""))
  pheatmap(mat, annotation_col = anno)
  dev.off()
  
  sampleDists <- dist(t(rld))
  sampleDists
  sampleDistMatrix <- as.matrix(sampleDists)
  pdf(paste("Sample_distance.rld",".(",condition,").","miR",".pdf", sep = ""))
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           annotation_col = anno)
  dev.off()

  res.df.05 <- as.matrix(res.df.05)
  rownames(res.df.05) <- sub(rownames(res.df.05), pattern = "\\..*$", replacement = "")
  rownames(res.df.05) <- mapIds(org.Hs.eg.db,
      keys=rownames(res.df.05),
      column="SYMBOL",
      keytype="ENSEMBL",
      multiVals="first")
  
res.df.05 <- res.df.05[!is.na(rownames(res.df.05)),]
res.df.05 <- as.data.frame(res.df.05[!duplicated(rownames(res.df.05)),])
res.df.05$symbol <- rownames(res.df.05)
geneList <- as.vector(res.df.05$padj)
names(geneList) <- as.vector(rownames(res.df.05))

sigGenes <- names(geneList)
sigmiRNAs <- rownames(data.miR.05)
sigmiRNAs.df <- as.data.frame(sigmiRNAs)

library(multiMiR)
multimir_res <- get_multimir(org = "hsa",
mirna = sigmiRNAs,
target = sigGenes,
table = "all",
summary = TRUE,
predicted.cutoff.type = "p",
predicted.cutoff = 20,
use.tibble = TRUE,
add.link = TRUE)
table(multimir_res@data$type)

for (val.group in c("validated", "predicted")) {
multimir_res.group <- paste('multimir_res', val.group, sep = ".")
multimir_res.tmp <- select(multimir_res, keytype = "type", keys = val.group, 
columns = columns(multimir_res))
assign(multimir_res.group, value = multimir_res.tmp)
rm(multimir_res.tmp)

multimir_res.df <- as.data.frame(get(multimir_res.group))
write.csv2(multimir_res.df,row.names=FALSE,file=paste("miRNA-mRNA_", val.group,"_interactions",".(",condition,")",".csv", sep = ""))
length(unique(get(multimir_res.group)$mature_mirna_id))
length(unique(get(multimir_res.group)$target_entrez))
unique.pairs <-
unique(data.frame(miRNA.ID = as.character(get(multimir_res.group)$mature_mirna_id),
target.symbol = as.character(get(multimir_res.group)$target_symbol)))
res.df.05$symbol <- toupper(rownames(res.df.05))
unique.pairs$target.symbol <- toupper(unique.pairs$target.symbol)
data.miR.05$symbol <- tolower(rownames(data.miR.05))
unique.pairs$miRNA.ID <- tolower(unique.pairs$miRNA.ID)
all.res.genes <- merge(x=unique.pairs, y=res.df.05, by.x="target.symbol", by.y="symbol", all.x=TRUE)
all.res.full <- merge(x=all.res.genes, y=data.miR.05, by.x="miRNA.ID", by.y="symbol", all.x=TRUE)
#all.res.full.p.001 <- subset(all.res.full,subset=padj<0.001 & QValue<0.001)

library(DESeq2)
mrna.dds.counts <- assays(dds)$counts
colnames(mrna.dds.counts) <- gsub(colnames(mrna.dds.counts),pattern="(i.*|_.*| ?ZBK|RymK ?)", replacement="")

edgeR.dge.mrna <- DGEList(counts=mrna.dds.counts)
#keep <- filterByExpr(y)
#y <- y[keep, , keep.lib.sizes=FALSE]
#normalization
#y <- estimateDisp(y, design)
#group <- factor(c(1,1,2,2,3,3))
#design <- model.matrix(~group)
#fit <- glmQLFit(y, design)
#qlf <- glmQLFTest(fit, coef=2:3)
#topTags(qlf)

miR.counts.norm.full.df <- as.data.frame(cpm(y, log = FALSE))
miR.counts.norm.full.df$group <- "miR"
miR.counts.norm.full.df.t <- t(miR.counts.norm.full.df)

mrna.dds.counts.norm.full.df <- as.data.frame(cpm(edgeR.dge.mrna, log = FALSE))
mrna.dds.counts.norm.full.df$group <- "mRNA"
mrna.dds.counts.norm.full.df.t <- t(mrna.dds.counts.norm.full.df)
merged.counts.norm.full.df.t <- merge(x = miR.counts.norm.full.df.t, y = mrna.dds.counts.norm.full.df.t, by.x = 0, by.y = 0)
rownames(merged.counts.norm.full.df.t) <- merged.counts.norm.full.df.t$Row.names
merged.counts.norm.full.df <- as.data.frame(t(merged.counts.norm.full.df.t))
merged.counts.norm.full.df <- merged.counts.norm.full.df[rownames(merged.counts.norm.full.df) != "Row.names",]
miR.counts.norm.full.df <- subset(merged.counts.norm.full.df, group == "miR")
mrna.dds.counts.norm.full.df <- subset(merged.counts.norm.full.df, group == "mRNA")

mrna.dds.counts.norm.full.df <- as.matrix(mrna.dds.counts.norm.full.df)
rownames(mrna.dds.counts.norm.full.df) <- sub(rownames(mrna.dds.counts.norm.full.df), pattern = "\\..*$", replacement = "")
rownames(mrna.dds.counts.norm.full.df) <- mapIds(org.Hs.eg.db,
       keys=rownames(mrna.dds.counts.norm.full.df),
       column="SYMBOL",
       keytype="ENSEMBL",
       multiVals="first")
mrna.dds.counts.norm.full.df <- mrna.dds.counts.norm.full.df[!is.na(rownames(mrna.dds.counts.norm.full.df)),]
mrna.dds.counts.norm.full.df <- mrna.dds.counts.norm.full.df[!duplicated(rownames(mrna.dds.counts.norm.full.df)),]

combined.df <- as.data.frame(t(merge(x=t(mrna.dds.counts.norm.full.df),y=t(miR.counts.norm.full.df), by.x="row.names", by.y="row.names")), stringsAsFactors=FALSE)
colnames(combined.df) <- unlist(combined.df["Row.names",])
combined.df <- combined.df[ -c(1),]
rownames(combined.df)[ combined.df$group == "miR" ] <- tolower(rownames(combined.df)[ combined.df$group == "miR" ])
rownames(combined.df)[ combined.df$group == "mRNA" ] <- toupper(rownames(combined.df)[ combined.df$group == "mRNA" ])

library(Hmisc)
correlation.df <- as.data.frame(cbind(c(""),c(""),c(""),c("")))
colnames(correlation.df) <- (c("mRNA.gene","miR.gene","R2","Pearson_p.value"))
correlation.df <- correlation.df[-c(1),]
for(i in seq(from=1, to=nrow(all.res.full))) {
	col1 <- all.res.full[i,c("target.symbol")]
	col2 <- all.res.full[i,c("miRNA.ID")]
	sub1 <- subset(combined.df,subset=rownames(combined.df) %in% col1 | rownames(combined.df) %in% col2)
	sub1$group <- NULL
	sub1 <- t(sub1)
	sub1 <- as.data.frame(sub1,stringsAsFactors=FALSE)
	sub1 <- sapply(sub1,as.numeric)
	if(ncol(sub1) != 2) {
		save.image("Error.RData")
		stop("The number of columns is not 2.")}
	cor1 <- rcorr(sub1)
	df1 <- as.data.frame(t(c(colnames(cor1$r),cor1$r[1,2], cor1$P[1,2])),stringsAsFactors=FALSE)
	colnames(df1) <-  colnames(correlation.df)
	correlation.df <- rbind(correlation.df,df1)}	
	correlation.df[,3] <- as.numeric(correlation.df[,3])
	correlation.df[,4] <- as.numeric(correlation.df[,4])
	all.res.full_correlation.df <- cbind(all.res.full, correlation.df)
	all.res.full_correlation.df <- all.res.full_correlation.df [order(all.res.full_correlation.df$Pearson_p.value),]

	
		write.csv2(all.res.full_correlation.df,row.names=FALSE,file=paste("miRNA-mRNA_unique_", val.group,"_interactions_with_Pearson_correlation_results",".(",condition,")",".csv", sep = ""))
	pdf(paste("Gene expression correlation results.",condition,".", val.group, "_interactions",".pdf", sep = ""))
for(i in seq(from=1, to=nrow(all.res.full_correlation.df))) {
	col1 <- all.res.full_correlation.df[i,c("mRNA.gene")]
	col2 <- all.res.full_correlation.df[i,c("miR.gene")]
	sub1 <- subset(combined.df,subset=rownames(combined.df) %in% col1 | rownames(combined.df) %in% col2)
	sub1$group <- NULL
	sub1 <- t(sub1)
	sub1 <- as.data.frame(sub1,stringsAsFactors=FALSE)
	sub1 <- sapply(sub1,as.numeric)
	if(ncol(sub1) != 2) {
		save.image("Error.RData")
		stop("The number of columns is not 2.")}
	cor1 <- rcorr(sub1)
	sub1 <- as.data.frame(sub1,stringsAsFactors=FALSE)
	title1 <- paste0("Gene expression correlation: ", col1, " vs ", col2 )
	label1 <- paste0("Pearson's correlation: R2 = ",signif(cor1$r[1,2],3),", p = ", signif(cor1$P[1,2],3))
	plot1 <- ggplot(sub1, aes(.data[[col1]],.data[[col2]])) + geom_point(size=3, col="red") + geom_smooth(method="lm", formula="y~x") + labs(title=title1, subtitle=label1) + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
	print(plot1)
}
	dev.off()
}

  library(reticulate)
  library(xlsx)
  analysis_type <- "ora"
  py_config()
  mieaa <- import(module = 'mieaa')
  mieaa_api <- mieaa$API()
  mieaa_api$invalidate()
  mieaa_categories <- names(unlist(mieaa_api$get_enrichment_categories(mirna_type = "mirna", species = "hsa")))
  categories.df <- as.data.frame(t(as.data.frame(mieaa_api$get_enrichment_categories(mirna_type = "mirna", species = "hsa"))))
  categories.df$Database <- rownames(categories.df)
  categories.df <- categories.df[,c(2,1)]
  categories.df <- setNames(categories.df,nm = c("Database", "Description"))
  xlsfile <- paste(condition, "mieaa_results", analysis_type, "xls",sep = ".")
  unlink(xlsfile)
  write.xlsx2(sigmiRNAs.df, sheetName="Significant miRNAs", row.names=FALSE,file=xlsfile, append = TRUE)
  write.xlsx2(categories.df, sheetName="Categories", row.names=FALSE,file=xlsfile, append = TRUE)
  
  if (length(sigmiRNAs) < 100) {map.width <- 20} else {map.width <- round(0.2 * length(sigmiRNAs))}
  
  pdf(paste(condition, "mieaa_results", analysis_type, "pdf",sep = "."), width = map.width, height = 20)
  for (i in mieaa_categories) {
    if (analysis_type == "ora") {
      cat(paste("Over-representation analysis - category:", i, "...",sep = " "))
      mieaa_api$run_ora(test_set = sigmiRNAs, categories = i, mirna_type = "mirna", species = 'hsa')
    } else if (analysis_type == "gsea") {
      cat(paste("Gene set enrichment analysis - category:", i, "...",sep = " "))
      mieaa_api$run_gsea(test_set = sigmiRNAs, categories = i, mirna_type = "mirna", species = 'hsa')
    }
    parameters <- as.data.frame(mieaa_api$get_enrichment_parameters())
    df_name <- paste(condition, analysis_type, i, "mieaa_parameters", sep = ".")
    assign(x = df_name, value = parameters)
    json <- mieaa_api$get_results(results_format = 'json')
    if (length(json) == 0) {
      cat(" Empty list found.\n")
      mieaa_api$invalidate()
      Sys.sleep(10)
      next }
    df <- data.frame(matrix(unlist(json), nrow=length(json), byrow=T),stringsAsFactors = FALSE)
    cols = c('category', 'subcategory', 'enrichment', 'p-value', 'p-adjusted', 'q-value', 'expected', 'observed', 'mirnas/precursors')
    colnames(df) <- cols
    df$`q-value` <- as.numeric(df$'q-value')
    df <- df[order(df[,6]),]
    df$rank <- -log10(df$`q-value`)
    df_name <- paste(condition, analysis_type, i, "mieaa_results", sep = ".")
    
    library(wordcloud)
    wordcloud(words = df$subcategory, freq = df$rank,max.words = 100, colors = rainbow(n = 100), scale = c(2,0.00002))
    title(main = paste(df_name, "Frequency: -log10 of adjusted p-values", sep = "\n"))
    
    list1 <- strsplit(as.matrix(df$`mirnas/precursors`),split = ";")
    
    if (length(list1) > 100) {
      list1 <- list1[1:100]
      list.100 = TRUE
      df.100 <- df[1:100,]
    }
    for (j in seq(length(list1))) {
      r1 <- rep(df$rank[j],length(list1[[j]]))
      names(r1) <- gsub(list1[[j]], pattern = " ", replacement = "")
      row_n <- paste0("row_", j)
      r1 <- as.data.frame(r1)
      r1 <- setNames(object = r1, nm = row_n)
      assign(x = row_n, value = r1)
      if (j == 1){
        all_rows <- get(row_n)
        next
      }
      all_rows <- merge(x = all_rows, y = get(row_n), by.x=0, by.y=0, all.x=TRUE, all.y=TRUE)
      rownames(all_rows) <- all_rows$Row.names
      all_rows$Row.names <- NULL
    }
    all_cols <- t(all_rows)
    if (exists('list.100')) {
      rownames(all_cols) <- df.100$subcategory } else {
        rownames(all_cols) <- df$subcategory }
    if (nrow(all_cols) > 2) {
      try(expr = pheatmap(all_cols, cluster_rows = FALSE, cluster_cols = FALSE, color = brewer.pal(n = 9, name = "Blues"), main = paste((df_name),"(-log10 of adjusted p-values)", sep = "\n")))
    }
    write.xlsx2(df, sheetName=i, row.names=FALSE,file=xlsfile, append = TRUE)
    assign(x = df_name, value = df)
    if (exists('list.100')) {
      rm(list.100) }
    cat(" Done.",sep="\n")
    mieaa_api$invalidate()
    Sys.sleep(10)
  }
  dev.off()

save.image("miR_mRNA.correlation.RData")

sessionInfo()
proc.time()

cat("The analysis is complete.\n")
