#! /usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
print("Program path: /shares/NGS/MichalinaZ_mRNA/gene_expression.mod.r")
if(length(arguments) != 3) {stop("The CSV file, exp. analysis type (DESEQ2, EDGER) and filtering status (FILTERED, UNFILTERED) were not provided.")}
csvfile <- arguments[1]
variant <- gsub(sub(csvfile, pattern = "\\.csv$", replacement = ""), pattern = ".*/", replacement = "")

args <- seq(1,14)
suffix <- arguments[2]
Filters <- arguments[3]

if(Filters == "FILTERED" & suffix == "DESEQ2") {
Filters <- "pOverA(0.20,log2(100));function(x)(IQR(x)>0.25)" }

origdir <- '/workspace/lukasz/Michalina_Z/New_analyses/ALL_BAMS'
workdir <- paste(origdir, variant, suffix, Filters, sep = "/")
dir.create(workdir, recursive = TRUE)
setwd(workdir)
getwd()
inputdir <- "/workspace/lukasz/Michalina_Z/New_analyses/ALL_BAMS/BAMS"
gtffile <- "/workspace/lukasz/NGS-all-in-one/REFSEQS/GRCh38.primary_assembly.annotation.gtf.gz"
threads <- 60
stranded <- "no"
alpha <- 0.1
lfcThreshold <- 0
sample_col <- "SampleID"
grouping_var <- "ALL_SAMPLES"
group_values <- "NA"
ind.factor <- "condition"
ind.factor2 <- NA

csvfile <- read.csv(csvfile, header=TRUE)
csvfile$SampleID <- gsub(gsub(gsub(csvfile$SampleID, pattern = "^[^0-9]+", replacement = ""), pattern = "[^0-9]+$", replacement = ""), pattern = "i?_[0-9]+", replacement = "")

filenames <- grep(list.files(path = inputdir, full.names = TRUE), pattern = ".*sorted.bam$", value = TRUE)
filenames.ost <- NULL
for(i in csvfile[[sample_col]]) {filenames.ost <- append(filenames.ost, grep(filenames, pattern = paste0("[^0-9]", i, "[^0-9]"), value = TRUE))}
filenames <- filenames.ost

#inputdir <- readLines(file("stdin"))

save.image(file=paste0(suffix,".RData"))
if(length(args)<14) {
	save.image(file=paste0(suffix,".RData"))
	stop("Fourteen arguments must be supplied in the following order (Expression dir, suffix, mapping dir, gtffile, strand specificity, alpha, LFC, filtering functions, threads, CSV file, sample column, grouping column, group, factor columns)", call.=FALSE) }

library("parallel")
library("BiocParallel")
register(MulticoreParam(workers=threads))
#register(SerialParam())

if (grouping_var != "ALL_SAMPLES") {
groups <- unlist(strsplit(grouping_var,split=","))
group.values <- unlist(strsplit(group_values, split=","))
num=0
for (group in groups)
{
num=num+1
group.column <- unlist(strsplit(group.values[grep(x=group.values, pattern=group)],split=":"))[1]
group.value <- unlist(strsplit(group.values[grep(x=group.values, pattern=group)],split=":"))[2]
name <- paste("df.subset",num, sep=".")
assign(x=name,value=subset(csvfile, csvfile[[group.column]]==group.value))
}
if (length(groups) == 2) {
	csvfile <- merge(x=df.subset.1, y=df.subset.2)
} else if (length(groups) == 1) {
	csvfile <- df.subset.1
} else { 
	save.image(file=paste0(suffix,".RData"))
	stop("The number of grouping variables cannot be higher than two.")}
}
if (length(args) == 14)
		{ csvfile <- subset(csvfile,!is.na(csvfile[[ind.factor]]))
	} else if (length(args) == 15)
		{ csvfile <- subset(csvfile,!is.na(csvfile[[ind.factor]]) & !is.na(csvfile[[ind.factor2]]))}
if (length(csvfile[[sample_col]])<2) {
	save.image(file=paste0(suffix,".RData"))
	stop("There are not enough samples to perform further analyses.")}

	library("Rsamtools")
	bamfiles <- BamFileList(filenames)
	library("GenomicFeatures")
	csvfile[[sample_col]] <- as.factor(csvfile[[sample_col]])
#	for(i in seq(from=14, to=length(args), by=1)){csvfile[[args[i]]] <- as.factor(csvfile[[args[i]]])}
	if(grepl(x = names(bamfiles[1]),pattern = "*.bam$")) {
	  bamsamples <- gsub(x=names(bamfiles),pattern=".sorted.bam$",replacement="")
	  bamsamples <- gsub(gsub(gsub(bamsamples, pattern = "^[^0-9]+", replacement = ""), pattern = "[^0-9]+$", replacement = ""), pattern = "i?_[0-9]+", replacement = "")
	  bamsamples <- as.data.frame(bamsamples)}
	names(bamsamples) <- c("names")
	sampleTable <- merge(x=csvfile, y=bamsamples, by.x=sample_col, by.y="names")
	rownames(sampleTable) <- sampleTable[[sample_col]]
		if (nrow(sampleTable)<2) {
		save.image(file=paste0(suffix,".RData"))
		stop("At least two samples are necessary to perform the gene expression analysis.")
		}
f.topGO <- function(geneList) {
  library("topGO")
  colMap <- function(x) {
    .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
    return(.col[match(1:length(x), order(x))])
  }
  GOdata.BP <- new("topGOdata",
                   ontology = "BP",
                   allGenes = geneList,
                   geneSel = topDiffGenes,
                   annotationFun=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
                   nodeSize = 10)
  if(length(GOdata.BP@graph@nodes) >= 50) {nodesn.BP <- 50} else {nodesn.BP <- length(GOdata.BP@graph@nodes)}
  
  resultFisher.BP <- runTest(GOdata.BP, algorithm = "classic", statistic = "fisher")
  resultFisher.BP
  resultKS.BP <- runTest(GOdata.BP, algorithm = "classic", statistic = "ks")
  resultKS.BP
  resultKS.elim.BP <- runTest(GOdata.BP, algorithm = "elim", statistic = "ks")
  resultKS.BP
  allRes.BP.classicKS <- GenTable(GOdata.BP, classicFisher = resultFisher.BP,
                                  classicKS = resultKS.BP, elimKS = resultKS.elim.BP,
                                  orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.BP)
  allRes.BP.elimKS <- GenTable(GOdata.BP, classicFisher = resultFisher.BP,
                               classicKS = resultKS.BP, elimKS = resultKS.elim.BP,
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes = nodesn.BP)
  allRes.BP.Fisher <- GenTable(GOdata.BP, classicFisher = resultFisher.BP,
                               classicKS = resultKS.BP, elimKS = resultKS.elim.BP,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.BP)
  geneSig <- geneList[topDiffGenes(geneList)]
  
  AnnotatedGenes.BP.classicKS <- sapply(allRes.BP.classicKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.BP, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.BP.classicKS[[x]][AnnotatedGenes.BP.classicKS[[x]] %in% names(geneSig)]}
  GeneList.BP.classicKS <- sapply(allRes.BP.classicKS$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.BP.classicKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.BP.classicKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.BP.classicKS)
  sink()
  
  AnnotatedGenes.BP.elimKS <- sapply(allRes.BP.elimKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.BP, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.BP.elimKS[[x]][AnnotatedGenes.BP.elimKS[[x]] %in% names(geneSig)]}
  GeneList.BP.elimKS <- sapply(allRes.BP.elimKS$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.BP.elimKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.BP.elimKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.BP.elimKS)
  sink()
  
  AnnotatedGenes.BP.Fisher <- sapply(allRes.BP.Fisher$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.BP, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.BP.Fisher[[x]][AnnotatedGenes.BP.Fisher[[x]] %in% names(geneSig)]}
  GeneList.BP.Fisher <- sapply(allRes.BP.Fisher$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.BP.Fisher",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.BP.Fisher",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.BP.Fisher)
  sink()
  
  GOdata.CC <- new("topGOdata",
                   ontology = "CC",
                   allGenes = geneList,
                   geneSel = topDiffGenes,
                   annotationFun=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
                   nodeSize = 10)
  if(length(GOdata.CC@graph@nodes) >= 50) {nodesn.CC <- 50} else {nodesn.CC <- length(GOdata.CC@graph@nodes)}
  
  resultFisher.CC <- runTest(GOdata.CC, algorithm = "classic", statistic = "fisher")
  resultFisher.CC
  resultKS.CC <- runTest(GOdata.CC, algorithm = "classic", statistic = "ks")
  resultKS.CC
  resultKS.elim.CC <- runTest(GOdata.CC, algorithm = "elim", statistic = "ks")
  resultKS.CC
  allRes.CC.classicKS <- GenTable(GOdata.CC, classicFisher = resultFisher.CC,
                                  classicKS = resultKS.CC, elimKS = resultKS.elim.CC,
                                  orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.CC)
  allRes.CC.elimKS <- GenTable(GOdata.CC, classicFisher = resultFisher.CC,
                               classicKS = resultKS.CC, elimKS = resultKS.elim.CC,
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes = nodesn.CC)
  allRes.CC.Fisher <- GenTable(GOdata.CC, classicFisher = resultFisher.CC,
                               classicKS = resultKS.CC, elimKS = resultKS.elim.CC,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.CC)
  geneSig <- geneList[topDiffGenes(geneList)]
  
  AnnotatedGenes.CC.classicKS <- sapply(allRes.CC.classicKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.CC, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.CC.classicKS[[x]][AnnotatedGenes.CC.classicKS[[x]] %in% names(geneSig)]}
  GeneList.CC.classicKS <- sapply(allRes.CC.classicKS$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.CC.classicKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.CC.classicKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.CC.classicKS)
  sink()
  
  AnnotatedGenes.CC.elimKS <- sapply(allRes.CC.elimKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.CC, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.CC.elimKS[[x]][AnnotatedGenes.CC.elimKS[[x]] %in% names(geneSig)]}
  GeneList.CC.elimKS <- sapply(allRes.CC.elimKS$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.CC.elimKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.CC.elimKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.CC.elimKS)
  sink()
  
  AnnotatedGenes.CC.Fisher <- sapply(allRes.CC.Fisher$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.CC, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.CC.Fisher[[x]][AnnotatedGenes.CC.Fisher[[x]] %in% names(geneSig)]}
  GeneList.CC.Fisher <- sapply(allRes.CC.Fisher$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.CC.Fisher",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.CC.Fisher",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.CC.Fisher)
  sink()
  
  GOdata.MF <- new("topGOdata",
                   ontology = "MF",
                   allGenes = geneList,
                   geneSel = topDiffGenes,
                   annotationFun=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
                   nodeSize = 10)
  if(length(GOdata.MF@graph@nodes) >= 50) {nodesn.MF <- 50} else {nodesn.MF <- length(GOdata.MF@graph@nodes)}
  
  resultFisher.MF <- runTest(GOdata.MF, algorithm = "classic", statistic = "fisher")
  resultFisher.MF
  resultKS.MF <- runTest(GOdata.MF, algorithm = "classic", statistic = "ks")
  resultKS.MF
  resultKS.elim.MF <- runTest(GOdata.MF, algorithm = "elim", statistic = "ks")
  resultKS.MF
  allRes.MF.classicKS <- GenTable(GOdata.MF, classicFisher = resultFisher.MF,
                                  classicKS = resultKS.MF, elimKS = resultKS.elim.MF,
                                  orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.MF)
  allRes.MF.elimKS <- GenTable(GOdata.MF, classicFisher = resultFisher.MF,
                               classicKS = resultKS.MF, elimKS = resultKS.elim.MF,
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes = nodesn.MF)
  allRes.MF.Fisher <- GenTable(GOdata.MF, classicFisher = resultFisher.MF,
                               classicKS = resultKS.MF, elimKS = resultKS.elim.MF,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.MF)
  geneSig <- geneList[topDiffGenes(geneList)]
  
  AnnotatedGenes.MF.classicKS <- sapply(allRes.MF.classicKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.MF, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.MF.classicKS[[x]][AnnotatedGenes.MF.classicKS[[x]] %in% names(geneSig)]}
  GeneList.MF.classicKS <- sapply(allRes.MF.classicKS$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.MF.classicKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.MF.classicKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.MF.classicKS)
  sink()
  
  AnnotatedGenes.MF.elimKS <- sapply(allRes.MF.elimKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.MF, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.MF.elimKS[[x]][AnnotatedGenes.MF.elimKS[[x]] %in% names(geneSig)]}
  GeneList.MF.elimKS <- sapply(allRes.MF.elimKS$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.MF.elimKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.MF.elimKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.MF.elimKS)
  sink()
  
  AnnotatedGenes.MF.Fisher <- sapply(allRes.MF.Fisher$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.MF, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.MF.Fisher[[x]][AnnotatedGenes.MF.Fisher[[x]] %in% names(geneSig)]}
  GeneList.MF.Fisher <- sapply(allRes.MF.Fisher$GO.ID,genes_with_GO_term)
  if (length(args)==14) {
    sink(paste("GO_top50-significant_genes.MF.Fisher",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.MF.Fisher",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.MF.Fisher)
  sink()
  
  library(xlsx)
  if (length(args)==14) {
    write.xlsx2(allRes.BP.classicKS,row.names=FALSE,file=paste("GO-top_50.BP",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "classicKS")
    write.xlsx2(allRes.BP.elimKS,row.names=FALSE,file=paste("GO-top_50.BP",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "elimKS", append = T)
    write.xlsx2(allRes.BP.Fisher,row.names=FALSE,file=paste("GO-top_50.BP",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "Fisher", append = T)
    
    write.xlsx2(allRes.CC.classicKS,row.names=FALSE,file=paste("GO-top_50.CC",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "classicKS")
    write.xlsx2(allRes.CC.elimKS,row.names=FALSE,file=paste("GO-top_50.CC",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "elimKS", append = T)
    write.xlsx2(allRes.CC.Fisher,row.names=FALSE,file=paste("GO-top_50.CC",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "Fisher", append = T)
    
    write.xlsx2(allRes.MF.classicKS,row.names=FALSE,file=paste("GO-top_50.MF",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "classicKS")
    write.xlsx2(allRes.MF.elimKS,row.names=FALSE,file=paste("GO-top_50.MF",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "elimKS", append = T)
    write.xlsx2(allRes.MF.Fisher,row.names=FALSE,file=paste("GO-top_50.MF",".(",ind.factor,").",suffix,".xls", sep = ""), sheetName = "Fisher", append = T)
  } else {
    write.xlsx2(allRes.BP.classicKS,row.names=FALSE,file=paste("GO-top_50.BP",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "classicKS")
    write.xlsx2(allRes.BP.elimKS,row.names=FALSE,file=paste("GO-top_50.BP",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "elimKS", append = T)
    write.xlsx2(allRes.BP.Fisher,row.names=FALSE,file=paste("GO-top_50.BP",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "Fisher", append = T)
    
    write.xlsx2(allRes.CC.classicKS,row.names=FALSE,file=paste("GO-top_50.CC",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "classicKS")
    write.xlsx2(allRes.CC.elimKS,row.names=FALSE,file=paste("GO-top_50.CC",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "elimKS", append = T)
    write.xlsx2(allRes.CC.Fisher,row.names=FALSE,file=paste("GO-top_50.CC",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "Fisher", append = T)
    
    write.xlsx2(allRes.MF.classicKS,row.names=FALSE,file=paste("GO-top_50.MF",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "classicKS")
    write.xlsx2(allRes.MF.elimKS,row.names=FALSE,file=paste("GO-top_50.MF",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "elimKS", append = T)
    write.xlsx2(allRes.MF.Fisher,row.names=FALSE,file=paste("GO-top_50.MF",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = ""), sheetName = "Fisher", append = T)
  }
  
  if (length(args)==14) {
    pdf(paste("GO-p-values_comparison",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-p-values_comparison",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  pValue.classic <- score(resultKS.BP)
  pValue.elim <- score(resultKS.elim.BP)[names(pValue.classic)]
  gstat <- termStat(GOdata.BP, names(pValue.classic))
  gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  gCol <- colMap(gstat$Significant)
  plot(x=pValue.classic, y=pValue.elim, xlab = "p-value - KS test (classic)", ylab = "p-value - KS test (elim)", pch = 19, cex = gSize, col = gCol, title(main="Gene ontology analysis (Biological process)"))
  
  pValue.classic <- score(resultKS.CC)
  pValue.elim <- score(resultKS.elim.CC)[names(pValue.classic)]
  gstat <- termStat(GOdata.CC, names(pValue.classic))
  gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  gCol <- colMap(gstat$Significant)
  plot(x=pValue.classic, y=pValue.elim, xlab = "p-value - KS test (classic)", ylab = "p-value - KS test (elim)", pch = 19, cex = gSize, col = gCol, title(main="Gene ontology analysis (Cellular component)"))
  
  pValue.classic <- score(resultKS.MF)
  pValue.elim <- score(resultKS.elim.MF)[names(pValue.classic)]
  gstat <- termStat(GOdata.MF, names(pValue.classic))
  gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  gCol <- colMap(gstat$Significant)
  plot(x=pValue.classic, y=pValue.elim, xlab = "p-value - KS test (classic)", ylab = "p-value - KS test (elim)", pch = 19, cex = gSize, col = gCol, title(main="Gene ontology analysis (Molecular function)"))
  dev.off()
  
  if (length(args)==14) {
    pdf(paste("GO-diagrams_top_10_test_KS_classic",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-diagrams_top_10_test_KS_classic",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  tryCatch(expr={showSigOfNodes(GOdata.BP, score(resultKS.BP), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (biological processes)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (biological processes): ERROR")}})	
  tryCatch(expr={showSigOfNodes(GOdata.CC, score(resultKS.CC), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (cellular components)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (cellular components): ERROR")}})
  tryCatch(expr={showSigOfNodes(GOdata.MF, score(resultKS.MF), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (molecular functions)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (molecular functions): ERROR")}})
  dev.off()
  
  if (length(args)==14) {
    pdf(paste("GO-diagrams_top_10_test_KS_elim",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-diagrams_top_10_test_KS_elim",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  tryCatch(expr={showSigOfNodes(GOdata.BP, score(resultKS.elim.BP), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (biological processes)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (biological processes): ERROR")}})	
  tryCatch(expr={showSigOfNodes(GOdata.CC, score(resultKS.elim.CC), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (cellular components)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (cellular components): ERROR")}})
  tryCatch(expr={showSigOfNodes(GOdata.MF, score(resultKS.elim.MF), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (molecular functions)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (molecular functions): ERROR")}})
  dev.off()
  
  if (length(args)==14) {
    pdf(paste("GO-diagrams_top_10_test_Fisher",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-diagrams_top_10_test_Fisher",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  tryCatch(expr={showSigOfNodes(GOdata.BP, score(resultFisher.BP), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (biological processes)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (biological processes): ERROR")}})	
  tryCatch(expr={showSigOfNodes(GOdata.CC, score(resultFisher.CC), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (cellular components)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (cellular components): ERROR")}})
  tryCatch(expr={showSigOfNodes(GOdata.MF, score(resultFisher.MF), firstSigNodes = 10, useInfo ='all'); title(main="Gene ontology (molecular functions)")},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="Gene ontology (molecular functions): ERROR")}})
  dev.off()
}
if(grepl(suffix, pattern="DESEQ2")) {
	if(! file.exists(paste0(suffix,".SE.transformed.filtered",".RData"))) {
		if (! file.exists(paste0(suffix,".SE.transformed",".RData"))) {

	if (length(args) == 14) { factors.formula <- as.formula(paste0("~",ind.factor)) } else { factors.formula <- as.formula(paste0("~",ind.factor,"+",ind.factor2)) }
	txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
	txdb
	tbg <- transcriptsBy(txdb, by="gene")
	save.image(file=paste0(suffix,".RData"))
	library("GenomicAlignments")
	if(stranded == "yes") {
		se <- summarizeOverlaps(features=tbg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=FALSE,
                        fragments=FALSE )
	} else {
		se <- summarizeOverlaps(features=tbg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE,
                        fragments=FALSE )	
	}
	se
	dim(se)
	assayNames(se)
	head(assay(se), 3)
	colSums(assay(se))
	rowRanges(se)
	str(metadata(rowRanges(se)))
	colnames(se)
	colData(se) <- DataFrame(sampleTable)
	colData(se)
	library("DESeq2")
	ds <- DESeqDataSet(se, design = factors.formula)
	if (length(ds[[sample_col]]) <= 30)
        {
	rld.e <- tryCatch(rld <- rlog(ds, blind = FALSE),error=function (e) { rlog(ds, blind = TRUE) })
	if (exists("rld")) { rm(rld.e) } else { rld <- rld.e}
	} else if (length(ds[[sample_col]]) > 30)
	{
	vsd <- vst(ds, blind = FALSE)
	}
	save.image(file=paste0(suffix,".SE.transformed",".RData"))
} else { load(file=paste0(suffix,".SE.transformed",".RData")) } 
library("DESeq2")
	if (exists("rld")) 
	{
		library(genefilter)
		if (Filters != "UNFILTERED")
		{
		ffunctions <- unlist(strsplit(Filters,split=";"))
		for(i in seq(from=1, to=length(ffunctions),by=1)) { name <- paste("ff",i,sep=".")
		assign(name,eval(parse(text=ffunctions[i])))}
		ff.number <- length(grep(objects(),pattern="ff.[0-9]+"))
			if (ff.number == 1) 
			{ selGenes <- genefilter(assay(rld),filterfun(ff.1))
		       	} else if (ff.number == 2) 
			{ selGenes <- genefilter(assay(rld),filterfun(ff.1,ff.2))
		       	} else if (ff.number == 3)
			{ selGenes <- genefilter(assay(rld),filterfun(ff.1,ff.2,ff.3))
		       	} else { 
				save.image(file=paste0(suffix,".SE.transformed",".RData"))
				stop("The number of filtering functions is incorrect.") }
			eset <- assay(rld[selGenes,])
			rld[rownames(assay(rld)) %in% rownames(eset), ] -> rld
			ds[rownames(assay(ds)) %in% rownames(eset), ] -> ds
			nrow(assay(ds))
		}
		save.image(file=paste0(suffix,".SE.transformed.filtered",".RData"))
	} else if(exists("vsd")) 
		{
		library(genefilter)
		if (Filters != "UNFILTERED")
		{
		ffunctions <- unlist(strsplit(Filters,split=";"))
		for(i in seq(from=1, to=length(ffunctions),by=1)) { name <- paste("ff",i,sep=".")
		assign(name,eval(parse(text=ffunctions[i])))}
		ff.number <- length(grep(objects(),pattern="ff.[0-9]+"))
			if (ff.number == 1) 
			{ selGenes <- genefilter(assay(vsd),filterfun(ff.1))
		       	} else if (ff.number == 2) 
			{ selGenes <- genefilter(assay(vsd),filterfun(ff.1,ff.2))
		       	} else if (ff.number == 3)
			{ selGenes <- genefilter(assay(vsd),filterfun(ff.1,ff.2,ff.3))
		       	} else { 
				save.image(file=paste0(suffix,".SE.transformed",".RData"))
				stop("The number of filtering functions is incorrect.") }
			eset <- assay(vsd[selGenes,])
			vsd[rownames(assay(vsd)) %in% rownames(eset), ] -> vsd
			ds[rownames(assay(ds)) %in% rownames(eset), ] -> ds
			nrow(assay(ds))
		}
		save.image(file=paste0(suffix,".SE.transformed.filtered",".RData"))
		}	
} else { load(file=paste0(suffix,".SE.transformed.filtered",".RData")) }
library(DESeq2)
RUNID <- gsub(x=suffix, pattern="(^[^.]*)(.*)", replacement="\\1")
for (i in grep(ls(), pattern="(^rld$|^vsd$)", value=TRUE)) {

pdf(paste("Read counts histogram-",suffix,".pdf", sep = ""))
hist(counts(ds), col="gray", main=paste0("DESeq2 - histogram of cumulative read counts in the ", RUNID, " run."))
dev.off()

library(ggplot2)
library(pals)
if (length(args)==14) {
	pdf(paste("PCA-plot.",i,".(",ind.factor,").",suffix,".pdf", sep = ""))
	plot.pca <- plotPCA(get(i), intgroup = ind.factor, returnData=TRUE)
	percentVar <- round(100 * attr(plot.pca, "percentVar"))
	if(NROW(plot.pca$name)<=25) {
		ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = name, shape = group)) +
        	geom_point(size =3) +
        	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		labs(color=sample_col, shape=ind.factor, title = "Principal component analysis (PCA) plot") +
		coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
		theme(plot.title = element_text(hjust = 0.5))
		print(ggplot.1)
		dev.off()
	} else 	{
		ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = group)) +
        	geom_point(size =3) +
        	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		labs(color=ind.factor, title = "Principal component analysis (PCA) plot") +
		coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
		theme(plot.title = element_text(hjust = 0.5))
		print(ggplot.1)
		dev.off()
		}
} else {
	pdf(paste("PCA-plot.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
	plot.pca <- plotPCA(get(i), intgroup = c(ind.factor,ind.factor2), returnData=TRUE)
	percentVar <- round(100 * attr(plot.pca, "percentVar"))
	if(NROW(plot.pca$name)<=25) {
		ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = name, shape = group)) +
        	geom_point(size =3) +
        	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		labs(color=sample_col, shape=paste(ind.factor,ind.factor2,sep=":"), title = "Principal component analysis (PCA) plot") +
		coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
		theme(plot.title = element_text(hjust = 0.5))
		print(ggplot.1)
		dev.off()
	} else 	{
		ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = group)) +
        	geom_point(size =3) +
        	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		labs(color=paste(ind.factor,ind.factor2,sep=":"), title = "Principal component analysis (PCA) plot") +
		coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
		theme(plot.title = element_text(hjust = 0.5))
		print(ggplot.1)
		dev.off()
		}
}

library("genefilter")
#library(tidyverse)
library("AnnotationDbi")
library("org.Hs.eg.db")

topVarGenes <- head(order(rowVars(assay(get(i))), decreasing = TRUE), 20)
mat  <- assay(get(i))[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
mat.ens <- mat 
rownames(mat) <- gsub(rownames(mat),pattern="\\..*$", replacement="")
rownames(mat) <- as.vector(mapIds(org.Hs.eg.db,
keys=rownames(mat),
column="SYMBOL",
keytype="ENSEMBL",
multiVals="first"))
for(j in seq(from=1, to=nrow(mat))) {rownames(mat)[j] [is.na(rownames(mat)[j])] <- rownames(mat.ens)[j]}

library("pheatmap")
library("RColorBrewer")

if (length(args)==14) {
	anno <- as.data.frame(colData(get(i))[, ind.factor],row.names = as.vector(sampleTable[[sample_col]]))
	colnames(anno) <- ind.factor
	pdf(paste("TopVarGenes.",i,".(",ind.factor,").",suffix,".pdf", sep = ""))
} else {
	anno <- as.data.frame(colData(get(i))[, c(ind.factor,ind.factor2)],row.names = as.vector(sampleTable[[sample_col]]))
        colnames(anno) <- c(ind.factor, ind.factor2)
	pdf(paste("TopVarGenes.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
}
pheatmap(mat, annotation_col = anno)
dev.off()

sampleDists <- dist(t(assay(get(i))))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )

#colnames(sampleDistMatrix) <- NULL
if (length(args)==14) {
  pdf(paste("Sample_distance.",i,".(",ind.factor,").",suffix,".pdf", sep = ""))
} else {
  pdf(paste("Sample_distance.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
}

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = anno
)
dev.off()
}
	library("DESeq2")
	if (length(ds[[sample_col]])<2) {
	save.image(file=paste0(suffix,".SE.transformed.filtered",".RData"))
    	stop('The analysis cannot be performed for less than two samples.')
    } 	else if (exists("rld.e"))
	{
    	res <- data.frame(
    	assay(rld),
    	avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
    	rLogFC = assay(rld)[,2] - assay(rld)[,1], FC = 2^(assay(rld)[,2] - assay(rld)[,1]) )
	res$Ensembl.ID <- gsub("\\..*","",rownames(res))
	library("AnnotationDbi")
	library("org.Hs.eg.db")
	res$symbol <- mapIds(org.Hs.eg.db,
	keys=res$Ensembl.ID,
	column="SYMBOL",
	keytype="ENSEMBL",
	multiVals="first")
	res$Entrez <- mapIds(org.Hs.eg.db,
	keys=res$Ensembl.ID,
	column="ENTREZID",
	keytype="ENSEMBL",
	multiVals="first")
	} else
	{
	dds <- DESeq(ds)
	pdf(paste("Read counts normalization-",suffix,".pdf", sep = ""))
	par(mfrow=c(2,1))
	boxplot(counts(dds)+1, col = "lightblue", las = 2, cex.names = 1, log="y")
	title(main=paste0("DESeq2 - raw read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
	boxplot(counts(dds, normalized=TRUE)+1, col = "lightblue", las = 2, cex.names = 1, log="y")
	title(main=paste0("DESeq2 - normalized read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
	dev.off()

    	res <- results(dds, alpha = alpha, lfcThreshold = lfcThreshold)
#	library(tidyverse)
	res$FC <- 2**(res$log2FoldChange)
	res$Ensembl.ID <- gsub("\\..*","",rownames(res))
	library("AnnotationDbi")
	library("org.Hs.eg.db")
	res$symbol <- mapIds(org.Hs.eg.db,
	keys=res$Ensembl.ID,
	column="SYMBOL",
	keytype="ENSEMBL",
	multiVals="first")
	res$Entrez <- mapIds(org.Hs.eg.db,
	keys=res$Ensembl.ID,
	column="ENTREZID",
	keytype="ENSEMBL",
	multiVals="first")
	geneList <- as.vector(res$padj)
	names(geneList) <- as.vector(res$symbol)
	geneList <- geneList[!is.na(geneList)]
	topDiffGenes <- function(pvalue) {return(pvalue < 0.05) }
	if(sum(topDiffGenes(geneList)) > 0) 
	{
	  f.topGO(geneList)
	}
	}
head(res)
#summary(res)
#mcols(res, use.names = TRUE)
if (! exists("rld.e")) {
library("apeglm")
if (length(args)==14) {
	pdf(paste("MA-plot",".(",ind.factor,").",suffix,".pdf", sep = ""))
} else {
	pdf(paste("MA-plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
}
plotMA(res, ylim = c(-5, 5))
dev.off()
}
if (length(args)==14) {
	pdf(paste("P-value-histogram",".(",ind.factor,").",suffix,".pdf", sep = ""))
} else {
	pdf(paste("P-value-histogram",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
}
if (! exists("rld.e")) {
hist(res$padj, col = "grey50", border = "white")
} else if (exists("rld.e"))
{
hist(res$FC, col = "grey50", border = "white")
}
dev.off()

if (! exists("rld.e")) {
resOrdered <- res[order(res$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
if (length(args)==14) {
	write.table(resOrderedDF, file = paste("DESeq2_analysis_results",".(",ind.factor,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
	library("ReportingTools")
	htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""),
        reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""))
} else {
	write.table(resOrderedDF, file = paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
        library("ReportingTools")
        htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""),
        reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""))
}
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url, browser="firefox")
} else if (exists("rld.e"))
{
resOrdered <- res[order(res$FC),]
head(resOrdered)
if (length(args)==14) {
        write.table(resOrdered, file = paste("DESeq2_analysis_results",".(",ind.factor,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
        library("ReportingTools")
        htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""),
        reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""))
} else {
        write.table(resOrdered, file = paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
        library("ReportingTools")
        htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""),
        reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""))
}
publish(resOrdered, htmlRep)
url <- finish(htmlRep)
browseURL(url, browser="firefox")
}

normalized.counts <- counts(dds, normalized=TRUE)
rownames(normalized.counts) <- gsub(rownames(normalized.counts),pattern="\\..*$", replacement="")
rownames(normalized.counts) <- as.vector(mapIds(org.Hs.eg.db,
                                  keys=rownames(normalized.counts),
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first"))
normalized.counts <- normalized.counts[!is.na(rownames(normalized.counts)),]
normalized.counts <- normalized.counts[order(rownames(normalized.counts)),]
normalized.counts <- cbind(normalized.counts, rownames(normalized.counts))
colnames(normalized.counts)[ncol(normalized.counts)] <- "Gene.name"
normalized.counts <- normalized.counts[,c(ncol(normalized.counts),1:ncol(normalized.counts)-1)]
normalized.counts <- as.data.frame(normalized.counts)

library(xlsx)
if (length(args)==14) {
  write.xlsx2(normalized.counts,row.names=FALSE,file=paste("Normalized_counts",".(",ind.factor,").",suffix,".xls", sep = ""))
} else {
  write.xlsx2(normalized.counts,row.names=FALSE,file=paste("Normalized_counts",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = "")) }

print("ColSums unfiltered:")
print(colSums(assay(se)))
print("ColSums filtered:")
print(colSums(counts(ds)))
final.RData.name <- paste0(suffix,".SE.transformed.filtered",".RData")
save.image(file = final.RData.name)
} else if(grepl(suffix, pattern="EDGER")) {
	if (! file.exists(paste0(suffix,".SE.transformed",".RData"))) {
		if (! file.exists(paste0(suffix,".SE.RData"))) {
	library(Rsubread)
	if(stranded == "yes") {
                se <- featureCounts(files=filenames,annot.ext=gtffile,
                        isGTFAnnotationFile=TRUE,GTF.featureType="transcript",GTF.attrType="gene_id",
                        isPairedEnd = FALSE,
			tmpDir=tempdir(),
                        nthreads=threads,
                        strandSpecific = 1 )
        } else {
                se <- featureCounts(files=filenames,annot.ext=gtffile,
                        isGTFAnnotationFile=TRUE,GTF.featureType="transcript",GTF.attrType="gene_id",
                        isPairedEnd = FALSE,
			tmpDir=tempdir(),
                        nthreads=threads,
                        strandSpecific = 0 )
        }
	save.image(file=paste0(suffix,".RData"))
	colnames(se$counts) <- bamsamples$names
	library(edgeR)
	edgeR.dge <- DGEList(counts=se$counts)
	edgeR.dge <- calcNormFactors(edgeR.dge)
	if (Filters != "UNFILTERED")
	{
		keep <- filterByExpr(edgeR.dge)
		edgeR.counts.filtered <- edgeR.dge$counts[keep,]
		edgeR.dge <- DGEList(counts=edgeR.counts.filtered)
		edgeR.dge <- calcNormFactors(edgeR.dge)
	}
	save.image(file=paste0(suffix,".SE.RData")) } else { load(file=paste0(suffix,".SE.RData")) }
	for(i in colnames(edgeR.dge$counts)) {temp <- edgeR.dge$counts[,i]/(edgeR.dge$samples$lib.size[rownames(edgeR.dge$samples) %in% i]*edgeR.dge$samples$norm.factors[rownames(edgeR.dge$samples) %in% i])*1e6; assign(value=temp,x=paste("col",i, sep="."))}
	colList <- grep(objects(),pattern="^col\\..+$",value=TRUE)
	edgeR.counts.norm <- sapply(mget(colList),cbind)
	colnames(edgeR.counts.norm) <- colnames(edgeR.dge$counts)
	rownames(edgeR.counts.norm) <- rownames(edgeR.dge$counts)
	edgeR.counts.norm.df <- as.data.frame(edgeR.counts.norm, stringsAsFactors=FALSE)
	rm(list=grep(objects(),pattern="^col\\..+$",value=TRUE))

	sampleTable <- merge(x=edgeR.dge$samples, y=csvfile, by.x="row.names", by.y=sample_col)
	colnames(sampleTable)[colnames(sampleTable) %in% "Row.names"] <- sample_col
	if (length(args)==14) { 
		design <- model.matrix(~sampleTable[[ind.factor]]) } else {
		design <- model.matrix(~sampleTable[[ind.factor]] + sampleTable[[ind.factor2]]) }
	edgeR.dge <- estimateDisp(edgeR.dge, design)
	rld <- cpm(edgeR.dge, log=TRUE) ## Shows normalized results in a log2 scale
		 #cpm(edgeR.dge, log=FALSE) ## Shows normalized results in a linear scale
	save.image(file=paste0(suffix,".SE.transformed",".RData"))
} else {load(file=paste0(suffix,".SE.transformed",".RData"))}
library(edgeR)

RUNID <- gsub(x=suffix, pattern="(^[^.]*)(.*)", replacement="\\1")
        pdf(paste("Read counts histogram-",suffix,".pdf", sep = ""))
	hist(edgeR.dge$counts, col="gray", main=paste0("edgeR - histogram of cumulative read counts in the ", RUNID, " run."))
	dev.off()
if (Filters != "UNFILTERED")
{
	pdf(paste("Read counts normalization-",suffix,".pdf", sep = ""))
	par(mfrow=c(2,1))
	boxplot(edgeR.dge$counts+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
	title(main=paste0("edgeR - raw read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
	boxplot(edgeR.counts.norm.df+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
	title(main=paste0("edgeR - normalized read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
	dev.off()
} else {
	pdf(paste("Read counts normalization-",suffix,".pdf", sep = ""))
        par(mfrow=c(2,1))
        boxplot(edgeR.dge$counts+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
        title(main=paste0("edgeR - raw read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
        boxplot(edgeR.counts.norm.df+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
        title(main=paste0("edgeR - normalized read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
        dev.off()
}

	library(pals)
	library(ggplot2)
	plot.mds <- plotMDS(edgeR.dge)
	dev.off()
	unlink('Rplots.pdf')
	plot.mds.df <- as.data.frame(plot.mds$cmdscale.out)
	if (length(args)==14) {
	  pdf(paste("MDS-plot",".(",ind.factor,").",suffix,".pdf", sep = ""))
	  if(nrow(plot.mds.df)<=25) {
	    ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = sampleTable[[sample_col]], shape = sampleTable[[ind.factor]])) + 
	      geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
	      xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
	      labs(title = "Multidimensional scaling (MDS) plot", shape = ind.factor, color = sample_col) + 
	      theme(plot.title = element_text(hjust = 0.5))
	    print(ggplot.1)
	  } else 	{
	    ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = sampleTable[[ind.factor]])) + 
	      geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
	      xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
	      labs(title = "Multidimensional scaling (MDS) plot", color = ind.factor) + 
	      theme(plot.title = element_text(hjust = 0.5))
	    print(ggplot.1)
	  }
	  dev.off()
	} else {
	  pdf(paste("MDS-plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
	  if(nrow(plot.mds.df)<=25) {
	    ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = sampleTable[[sample_col]], shape = interaction(sampleTable[[ind.factor]], sampleTable[[ind.factor2]]))) + 
	      geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
	      xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
	      labs(title = "Multidimensional scaling (MDS) plot", shape = paste(ind.factor, ind.factor2, sep = "+"), color = sample_col) + 
	      theme(plot.title = element_text(hjust = 0.5))
	    print(ggplot.1)
	  } else 	{
	    ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = interaction(sampleTable[[ind.factor]], sampleTable[[ind.factor2]]))) + 
	      geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
	      xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
	      labs(title = "Multidimensional scaling (MDS) plot", color = paste(ind.factor, ind.factor2, sep = "+")) + 
	      theme(plot.title = element_text(hjust = 0.5))
	    print(ggplot.1)
	  }
	  dev.off()
	}
	
library("genefilter")
#library(tidyverse)
library("AnnotationDbi")
library("org.Hs.eg.db")

topVarGenes <- head(order(rowVars(rld), decreasing = TRUE), 20)
mat  <- rld[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
mat.ens <- mat
rownames(mat) <- gsub(rownames(mat),pattern="\\..*$", replacement="")
rownames(mat) <- as.vector(mapIds(org.Hs.eg.db,
keys=rownames(mat),
column="SYMBOL",
keytype="ENSEMBL",
multiVals="first"))
for(i in seq(from=1, to=nrow(mat))) {rownames(mat)[i] [is.na(rownames(mat)[i])] <- rownames(mat.ens)[i]}

library("pheatmap")
library("RColorBrewer")

if (length(args)==14) {
        anno <- as.data.frame(sampleTable[, ind.factor],row.names = as.vector(sampleTable[[sample_col]]))
        colnames(anno) <- ind.factor
        pdf(paste("TopVarGenes.rld",".(",ind.factor,").",suffix,".pdf", sep = ""))
} else {
        anno <- as.data.frame(sampleTable[, c(ind.factor, ind.factor2)],row.names = as.vector(sampleTable[[sample_col]]))
        colnames(anno) <- c(ind.factor, ind.factor2)
        pdf(paste("TopVarGenes.rld",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
}
pheatmap(mat, annotation_col = anno)
dev.off()

sampleDists <- dist(t(rld))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)

#colnames(sampleDistMatrix) <- NULL
if (length(args)==14) {
  pdf(paste("Sample_distance.rld",".(",ind.factor,").",suffix,".pdf", sep = ""))
} else {
  pdf(paste("Sample_distance.rld",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
}

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = anno
)
dev.off()

if (length(sampleTable[[sample_col]])<3) {
	save.image(file=paste0(suffix,".SE.transformed",".RData"))
    	stop('The analysis cannot be performed for less than three samples.')}
	fit <- glmQLFit(edgeR.dge, design)
	qlf <- glmQLFTest(fit)
	res <- topTags(qlf, n=nrow(qlf$table))
	res$table$FC <- 2**res$table$logFC
       	res$table$Ensembl.ID <- gsub("\\..*","",rownames(res))
	library("AnnotationDbi")
	library("org.Hs.eg.db")
	res$table$symbol <- mapIds(org.Hs.eg.db,
	keys=res$table$Ensembl.ID,
	column="SYMBOL",
	keytype="ENSEMBL",
	multiVals="first")
	res$table$Entrez <- mapIds(org.Hs.eg.db,
	keys=res$table$Ensembl.ID,
	column="ENTREZID",
	keytype="ENSEMBL",
	multiVals="first")
	geneList <- as.vector(res$table$FDR)
	names(geneList) <- as.vector(res$table$symbol)
	geneList <- geneList[!is.na(geneList)]
	topDiffGenes <- function(pvalue) {return(pvalue < 0.05) }
	if(sum(topDiffGenes(geneList)) > 0) 
	{
f.topGO(geneList)
	}
head(res$table)
	if (length(args)==14) {
		pdf(paste("P-value-histogram",".(",ind.factor,").",suffix,".pdf", sep = ""))
	} else {
		pdf(paste("P-value-histogram",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
	}
	hist(res$table$FDR, col = "grey50", border = "white")
	dev.off()

	resOrderedDF <- res$table[order(res$table$FDR),]
	head(resOrderedDF)
	if (length(args)==14) {
        	write.table(resOrderedDF, file = paste("edgeR_analysis_results",".(",ind.factor,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
        	library("ReportingTools")
        	htmlRep <- HTMLReport(shortName="edgeR analysis report", title=paste("edgeR_analysis_results",".(",ind.factor,").",suffix, sep = ""),
        	reportDirectory=paste("edgeR_analysis_results",".(",ind.factor,").",suffix, sep = ""))
	} else {
        	write.table(resOrderedDF, file = paste("edgeR_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
        	library("ReportingTools")
        	htmlRep <- HTMLReport(shortName="edgeR analysis report", title=paste("edgeR_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""),
        	reportDirectory=paste("edgeR_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""))
	}
	publish(resOrderedDF, htmlRep)
	url <- finish(htmlRep)
	browseURL(url, browser="firefox")
	
	normalized.counts <- cpm(edgeR.dge, log=FALSE)
	rownames(normalized.counts) <- gsub(rownames(normalized.counts),pattern="\\..*$", replacement="")
	rownames(normalized.counts) <- as.vector(mapIds(org.Hs.eg.db,
	                                                keys=rownames(normalized.counts),
	                                                column="SYMBOL",
	                                                keytype="ENSEMBL",
	                                                multiVals="first"))
	normalized.counts <- normalized.counts[!is.na(rownames(normalized.counts)),]
	normalized.counts <- normalized.counts[order(rownames(normalized.counts)),]
	normalized.counts <- cbind(normalized.counts, rownames(normalized.counts))
	colnames(normalized.counts)[ncol(normalized.counts)] <- "Gene.name"
	normalized.counts <- normalized.counts[,c(ncol(normalized.counts),1:ncol(normalized.counts)-1)]
	normalized.counts <- as.data.frame(normalized.counts)
	
	library(xlsx)
	if (length(args)==14) {
	  write.xlsx2(normalized.counts,row.names=FALSE,file=paste("Normalized_counts",".(",ind.factor,").",suffix,".xls", sep = ""))
	} else {
	  write.xlsx2(normalized.counts,row.names=FALSE,file=paste("Normalized_counts",".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = "")) }

	print("ColSums unfiltered:")
	print(colSums(se$counts))
	print("ColSums filtered:")
	print(colSums(edgeR.dge$counts))
	final.RData.name <- paste0(suffix,".SE.transformed",".RData")
	
	save.image(file = final.RData.name)
}
library("ReactomePA")
library("ggplot2")
res.df <- as.data.frame(res)
if (suffix == "DESEQ2") {
  res.df.05 <- subset(res.df, padj < 0.05) } else if (suffix == "EDGER") {
    res.df.05 <- subset(res.df, FDR < 0.05)}

for (exp.type in c("upregulated genes", "downregulated genes")) {
  if (exp.type == "upregulated genes") {
    res.df.05.FC1.5 <- subset(res.df.05, FC >= 3/2)
  } else {
    res.df.05.FC1.5 <- subset(res.df.05, FC <= 2/3)
  }
  geneList <- res.df.05.FC1.5$FC
  names(geneList) <- res.df.05.FC1.5$Entrez
  geneList <- geneList[!is.na(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- geneList[!duplicated(names(geneList))]
  geneList.df <- as.data.frame(geneList)
  colnames(geneList.df) <- "FC"
  if (nrow(geneList.df) > 0) {
    rownames(geneList.df) <- as.vector(mapIds(org.Hs.eg.db,
                                              keys=rownames(geneList.df),
                                              column="SYMBOL",
                                              keytype="ENTREZID",
                                              multiVals="first"))}
  geneList.df.name <- paste("geneList",exp.type, sep = ".")
  assign(geneList.df.name, value = geneList.df)
  
  de <- names(geneList)
  if (length(de) > 0) {
    de.name <- paste("de",exp.type, sep = ".")
    assign(de.name, value = de)
    
    temp.x <- enrichPathway(gene=get(de.name), organism = "human", pAdjustMethod = "BH", pvalueCutoff=0.05, readable=TRUE)
    x.name <- paste("x",exp.type, sep = ".")
    x.name.df <- paste("x",exp.type, "df", sep = ".")
    assign(x.name, value = temp.x)
    rm(temp.x)
    try(p1 <- dotplot(get(x.name), showCategory=50) + labs(title = paste("Pathway enrichment analysis", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    p1.name <- paste('p1', exp.type, sep = ".")
    if(nrow(p1$data) > 0) {assign(p1.name, value = p1)
      p1.pathways <- as.character(p1$data[order(p1$data$GeneRatio, decreasing = T),][["Description"]][1:if(nrow(p1$data)<3){nrow(p1$data)} else {3}])
      if(exp.type == "upregulated genes") {
      pdf(paste("Reactome analysis", suffix, "01_1", "pdf", sep = "."), height = 10, width = 13)} else if(exp.type == "downregulated genes") {
      pdf(paste("Reactome analysis", suffix, "07_1", "pdf", sep = "."), height = 10, width = 13)}
      for(i in p1.pathways) {try(expr = {pp1 <- viewPathway(i)
      pp1 <- pp1 + labs(title = paste(exp.type, i, sep = ": ")) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
      print(pp1)})}
      dev.off()
      }
    try(p2 <- heatplot(get(x.name), showCategory = 50, foldChange = geneList) + labs(title = paste("Pathway enrichment heatmap", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x  = element_line(color = "grey", linetype = "solid", size = 0.5), panel.grid.major.y  = element_line(color = "grey", linetype = "solid", size = 0.5)))
    p2.name <- paste('p2', exp.type, sep = ".")
    if(exists('p2')) {assign(p2.name, value = p2)}
    try(p3 <- emapplot(get(x.name), showCategory = 50, color = "p.adjust") + labs(title = paste("Pathway enrichment map", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    p3.name <- paste('p3', exp.type, sep = ".")
    if(exists('p3')) {assign(p3.name, value = p3)}
    try(p4 <- cnetplot(get(x.name), categorySize="qvalue", foldChange=geneList) + labs(title = paste("Complex associations map", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    p4.name <- paste('p4', exp.type, sep = ".")
    if(exists('p4')) {assign(p4.name, value = p4)}
    
    tmp.df <- as.data.frame(get(x.name))
    assign(x.name.df, value = tmp.df)
    
    try({y <- gsePathway(geneList, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH", by = "fgsea"); res.fgsea <- as.data.frame(y); p5 <- emapplot(y, showCategory = 50, color = "p.adjust") + labs(title = paste("Gene set enrichment analysis", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))})
    if(exists("y")) { if(nrow(y) > 0) {
      for(i in seq(1,length(y@result$core_enrichment))) {y@result$core_enrichment[i] <- paste(as.character(mapIds(org.Hs.eg.db, keys = unlist(strsplit(y@result$core_enrichment[i], split = "/")), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")),collapse = "/")}}}
    p5.name <- paste('p5', exp.type, sep = ".")
    if(exists('p5')) {assign(p5.name, value = p5)}
    
    y.name <- paste("y",exp.type, sep = ".")
    y.name.df <- paste("y",exp.type, "df", sep = ".")
    if(exists("y")) {assign(y.name, value = y)}
    if(exists("y")) {assign(y.name.df, value = as.data.frame(get(y.name)))}
    geneList.Symbols <- geneList
    if (length(geneList.Symbols) > 0) {
      names(geneList.Symbols) <- mapIds(org.Hs.eg.db, keys = names(geneList.Symbols), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")}
    if(exists("y")) {if(nrow(get(y.name.df)) > 0) { 
      p6 <- heatplot(get(y.name), showCategory = 50, foldChange = geneList.Symbols) + labs(title = paste("Gene set enrichment heatmap", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x  = element_line(color = "grey", linetype = "solid", size = 0.5), panel.grid.major.y  = element_line(color = "grey", linetype = "solid", size = 0.5))}}
    p6.name <- paste('p6', exp.type, sep = ".")
    if(exists('p6')) {assign(p6.name, value = p6)}
    rm(p1, p2, p3, p4, p5, p6, y)
  }}
if(! exists("de.upregulated genes")) {`de.upregulated genes` <- NULL}
if (! exists("de.downregulated genes")) {`de.downregulated genes` <- NULL}
de.list.full <- list(`de.upregulated genes`, `de.downregulated genes`)
names(de.list.full) <- c("upregulated genes", "downregulated genes")
library(clusterProfiler)
try(compareClustersRes <- compareCluster(de.list.full, fun="enrichPathway", organism = "human", pAdjustMethod = "BH", pvalueCutoff=0.05, readable=TRUE))
try(p7 <- dotplot(compareClustersRes, showCategory=50) + labs(title = paste("Pathway enrichment analysis - group comparison", sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
try(compareClustersRes.df <- as.data.frame(compareClustersRes))

if (exists("p1.upregulated genes")){
  pdf(paste("Reactome analysis", suffix, "01", "pdf", sep = "."), height = 10, width = 13)
  print(`p1.upregulated genes`)
  dev.off()}
if(exists('p2.upregulated genes')) { if(nrow(`p2.upregulated genes`$data) > 260){
  pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = nrow(`p2.upregulated genes`$data)/20)} else {
    pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = 13)
  }} else {
    pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = 10)}
if(exists('p2.upregulated genes')){
  print(`p2.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment heatmap: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
dev.off()

pdf(paste("Reactome analysis", suffix, "03", "pdf", sep = "."), height = 13, width = 13)
if(exists('p3.upregulated genes')) {
  print(`p3.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment map: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
dev.off()
pdf(paste("Reactome analysis", suffix, "04", "pdf", sep = "."), height = 10, width = 10)
if(exists('p4.upregulated genes')) {
  print(`p4.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Complex assocations map: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
dev.off()
pdf(paste("Reactome analysis", suffix, "05", "pdf", sep = "."), height = 10, width = 10)
if(exists('p5.upregulated genes')) {
  print(`p5.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment analysis: no term enriched under specific pvalueCutoff", "upregulated genes", sep = "-"))}
dev.off()
if(exists('p6.upregulated genes')) { if (nrow(`p6.upregulated genes`$data) > 260){
  pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = nrow(`p6.upregulated genes`$data)/20)} else {
    pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = 13)
  }} else {
    pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = 10)}
if(exists('p6.upregulated genes')) {
  print(`p6.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment heatmap: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
dev.off()
if (exists("p1.downregulated genes")){
  pdf(paste("Reactome analysis", suffix, "07", "pdf", sep = "."), height = 10, width = 13)
  print(`p1.downregulated genes`)
  dev.off()}
if(exists('p2.downregulated genes')) { if(nrow(`p2.downregulated genes`$data) > 260){
  pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = nrow(`p2.downregulated genes`$data)/20)} else {
    pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = 13)
  }} else {
    pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = 10)}
if(exists('p2.downregulated genes')){
  print(`p2.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment heatmap: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
dev.off()
pdf(paste("Reactome analysis", suffix, "09", "pdf", sep = "."), height = 13, width = 13)
if(exists('p3.downregulated genes')) {
  print(`p3.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment map: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
dev.off()
pdf(paste("Reactome analysis", suffix, "10", "pdf", sep = "."), height = 10, width = 10)
if(exists('p4.downregulated genes')) {
  print(`p4.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Complex assocations map: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
dev.off()
pdf(paste("Reactome analysis", suffix, "11", "pdf", sep = "."), height = 10, width = 10)
if(exists('p5.downregulated genes')) {
  print(`p5.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment analysis: no term enriched under specific pvalueCutoff", "downregulated genes", sep = "-"))}
dev.off()
if(exists('p6.downregulated genes')) { if(nrow(`p6.downregulated genes`$data) > 260){
  pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = nrow(`p6.downregulated genes`$data)/20)} else {
    pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width =  13)
  }} else {
    pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = 10)}
if(exists('p6.downregulated genes')){
  print(`p6.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment heatmap: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
dev.off()
if(exists("p7")) {
  if(nrow(p7$data)>60) {p7.height = floor(nrow(p7$data)/6)} else {p7.height = 10}
  pdf(paste("Reactome analysis", suffix, "13", "pdf", sep = "."), height = p7.height, width = 13)
  print(p7)
  dev.off()}

pdffiles <- sort(list.files(pattern = paste("Reactome analysis", suffix, ".*", "pdf", sep = ".")), method = "radix")
library(pdftools)
if (length(args)==14) {
  reactome.pdf.name = paste0("Reactome analysis", ".(",ind.factor,").",suffix,".pdf", sep = "")
} else {
  reactome.pdf.name = paste0("Reactome analysis", ".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = "")
}
pdf_combine(pdffiles, output = reactome.pdf.name)
unlink(pdffiles)

if (length(args)==14) {
  reactome.xls.name = paste0("Reactome analysis", ".(",ind.factor,").",suffix,".xls", sep = "")
} else {
  reactome.xls.name = paste0("Reactome analysis", ".(",ind.factor,"+",ind.factor2,").",suffix,".xls", sep = "")
}

unlink(reactome.xls.name)
if(exists("geneList.upregulated genes")) {write.xlsx2(`geneList.upregulated genes`, sheetName = "Upregulated.genes", file = reactome.xls.name, row.names = TRUE)}
if(exists("x.upregulated genes.df")) {write.xlsx2(`x.upregulated genes.df`, sheetName = "Path.enrich.upregulated.genes", file = reactome.xls.name, append = TRUE, row.names = FALSE)}
if(exists("y.upregulated genes.df")) {write.xlsx2(`y.upregulated genes.df`, sheetName = "Gene.set.enrich.upregulated.genes", file = reactome.xls.name, append = TRUE, row.names = FALSE)}
if(exists("geneList.downregulated genes")) {write.xlsx2(`geneList.downregulated genes`, sheetName = "Downregulated.genes", file = reactome.xls.name, append = TRUE, row.names = TRUE)}
if(exists("x.downregulated genes.df")) {write.xlsx2(`x.downregulated genes.df`, sheetName = "Path.enrich.downregulated.genes", file = reactome.xls.name, append = TRUE, row.names = FALSE)}
if(exists("y.downregulated genes.df")) {write.xlsx2(`y.downregulated genes.df`, sheetName = "Gene.set.enrich.downregulated.genes", file = reactome.xls.name, append = TRUE, row.names = FALSE)}
if(exists("compareClustersRes.df")) {write.xlsx2(compareClustersRes.df, sheetName = "Path.enrich.group.comparison", file = reactome.xls.name, append = TRUE, row.names = FALSE)}

save.image(file = final.RData.name)

cat("All done.\n")
