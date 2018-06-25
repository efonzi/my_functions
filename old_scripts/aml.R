
amlDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml"
circosTutorialDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/tutorial_papers/circos-tutorials-0.67/data"

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------- 170426 - CHROMOTHRIPSIS - ADD CHR INFO TO GENE FISHER AND CLUSTERING OF QVALS -----------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
library(org.Hs.eg.db)
library(ggplot2)

# annotate gene coordinates for table containing gene symbols
# "gene_col" is the number of the column containing the gene symbols
getCoordinatesForSymbols <- function(genes, gene_col){
  
  # get gene list
  symbol = genes[,gene_col]
  
  # get symbol2entrezID table
  genes2entrez = select(org.Hs.eg.db, keys=as.character(symbol), columns="ENTREZID", keytype = "SYMBOL")
  colnames(genes2entrez)[2] = "entrez"
  
  # as it has multiple entrezID per symbol, get those symbols
  sel = genes2entrez[which(duplicated(genes2entrez$SYMBOL)), "SYMBOL"]
  dup = genes[genes[,gene_col] %in% sel,]
  # then add rows for duplicated symbols to "genes"
  genes = rbind(genes, dup)
  
  # sort both tables by symbol and merge
  genes = genes[order(genes[,gene_col]),]
  genes2entrez = genes2entrez[order(genes2entrez$SYMBOL),]
  genes = cbind(genes, entrez=genes2entrez$entrez, stringsAsFactors=F)
  
  # remove symbols with unknown entrezID
  genes = genes[complete.cases(genes),]
  
  # add columns
  genes$chr = rep(NA, dim(genes)[1])
  genes$start = rep(NA, dim(genes)[1])
  genes$end = rep(NA, dim(genes)[1])
  genes$strand = rep("+", dim(genes)[1])
  genes$chrArm = rep(NA, dim(genes)[1])
  genes$cytoband = rep(NA, dim(genes)[1])
  
  # get list of start and end positions of entrez ids
  startPos = mget(genes$entrez,envir=org.Hs.egCHRLOC, ifnotfound=NA)
  endPos = mget(genes$entrez,envir=org.Hs.egCHRLOCEND, ifnotfound=NA)
  # get list of cytobands of entrez ids
  cytoband = mget(genes$entrez,envir=org.Hs.egMAP, ifnotfound=NA)
  
  # loop over table rows and fill start positions, end positions and strand sign of genes
  for(r in 1:dim(genes)[1]){
    if(r%%1000==0){cat("filling info for gene", r, "/", dim(genes)[1], "\n")}
    entrez = genes$entrez[r]
    if(!is.na(entrez)){
      if(!is.na(startPos[[entrez]][1])){
        genes$chr[r] = as.numeric(names(startPos[[entrez]][1]))
        genes$start[r] = as.character(startPos[[entrez]][1])
        genes$end[r] = as.character(endPos[[entrez]][1])
        genes$chrArm[r] = cytoband[[entrez]]
        genes$cytoband[r] = cytoband[[entrez]]
        if(grepl("-", genes$start[r])){
          genes$start[r] = as.numeric(strsplit(genes$start[r], "-")[[1]][2]) # positions on negative strands come as negative integers
          genes$end[r] = as.numeric(strsplit(genes$end[r], "-")[[1]][2]) # positions on negative strands come as negative integers
          genes$strand[r] = "-"}
      }
    }
  }
  
  # extract chromosome arm info
  genes[grep("q", genes$chrArm), "chrArm"] = as.vector(sapply(genes[grep("q", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "q")[[1]][1], "q")))
  genes[grep("p", genes$chrArm), "chrArm"] = as.vector(sapply(genes[grep("p", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "p")[[1]][1], "p")))
  
  # remove symbols with unknown chromosome
  genes = genes[!is.na(genes$chrArm),]
  
  genes$start = as.numeric(genes$start)
  genes$end = as.numeric(genes$end)
  
  return(genes)
  
}

# load sample list and select case and ctrl
sampleID = read.delim(paste0(amlDIR, "/chromothripsis/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
selCase = sampleID[sampleID$chromothripsis==1,] # get case samples
selCtrl = sampleID[sampleID$chromothripsis==2,] # get ctrl samples

# load results of fisher on genes for loss/del
lossDel = read.delim(paste0(amlDIR, "/chromothripsis/170414_geneFisher_lossDel.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)

# add entrezID and chromosome coordinates
lossDel = getCoordinatesForSymbols(lossDel, 1)

lossDel$cytoband2 = sapply(lossDel$cytoband, function(x) strsplit(x, "[.]")[[1]][1])
lossDel$cytoband2 = sapply(lossDel$cytoband2, function(x) strsplit(x, "-")[[1]][1])
lossDel$cytoband2 = sapply(lossDel$cytoband2, function(x) strsplit(x, " ")[[1]][1])

boxplot(lossDel$QVAL~lossDel$chrArm)

ld4 = lossDel[lossDel$QVAL < 10^-4,]
as.data.frame(sort(table(ld4$cytoband)))
as.data.frame(sort(table(ld4$cytoband2)))
ggplot(ld4, aes(QVAL, fill=chrArm)) + geom_histogram(bins=50)
ggplot(ld4, aes(log10(QVAL), fill=chrArm)) + geom_histogram(binwidth = 0.05)

ld6 = lossDel[lossDel$QVAL < 10^-6,]
ggplot(ld6, aes(log10(QVAL), fill=cytoband)) + geom_histogram(binwidth = 0.05)
ggplot(ld6, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05) + scale_fill_discrete(name="cytoband")

ld10 = lossDel[lossDel$QVAL < 10^-10,]
ggplot(ld10, aes(log10(QVAL), fill=cytoband)) + geom_histogram(binwidth = 0.05)
ggplot(ld10, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05)

ld14 = lossDel[lossDel$QVAL < 10^-14,]
ggplot(ld14, aes(log10(QVAL), fill=cytoband)) + geom_histogram(binwidth = 0.05)
ggplot(ld14, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05)

ld4_6 = lossDel[lossDel$QVAL < 10^-4 & lossDel$QVAL >= 10^-6,]
sort(table(ld4_6$cytoband))
sort(table(ld4_6$cytoband2))

ggplot(ld4_6, aes(log10(QVAL), fill=chrArm)) + geom_histogram(binwidth = 0.05)
ggplot(ld4_6, aes(log10(QVAL), fill=cytoband)) + geom_histogram(binwidth = 0.05)
ggplot(ld4_6, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05)

boxplot(ld4$QVAL~ld4$chrArm)
boxplot(ld4$QVAL~ld4$cytoband2)
boxplot(ld10$QVAL~ld10$cytoband2)


write.table(lossDel, file=paste0(amlDIR, "/chromothripsis/170426_geneFisher_lossDel_chr_cytoband.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)

jpeg(filename=paste0(amlDIR, "/chromothripsis/170426_lossDel_qval_histogram_4.jpeg"), width=9, height=6, units="in", res=175, quality=100)
ggplot(ld4, aes(log10(QVAL), fill=chrArm)) + geom_histogram(binwidth = 0.05)
dev.off()
jpeg(filename=paste0(amlDIR, "/chromothripsis/170426_lossDel_qval_histogram_6.jpeg"), width=9, height=6, units="in", res=175, quality=100)
ggplot(ld6, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05) + scale_fill_discrete(name="cytoband")
dev.off()
jpeg(filename=paste0(amlDIR, "/chromothripsis/170426_lossDel_qval_histogram_10.jpeg"), width=9, height=6, units="in", res=175, quality=100)
ggplot(ld10, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05) + scale_fill_discrete(name="cytoband")
dev.off()
jpeg(filename=paste0(amlDIR, "/chromothripsis/170426_lossDel_qval_histogram_14.jpeg"), width=9, height=6, units="in", res=175, quality=100)
ggplot(ld14, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05) + scale_fill_discrete(name="cytoband")
dev.off()
jpeg(filename=paste0(amlDIR, "/chromothripsis/170426_lossDel_qval_histogram_4-6.jpeg"), width=9, height=6, units="in", res=175, quality=100)
ggplot(ld4_6, aes(log10(QVAL), fill=chrArm)) + geom_histogram(binwidth = 0.05)
dev.off()
jpeg(filename=paste0(amlDIR, "/chromothripsis/170426_lossDel_qval_histogram_4-6.jpeg"), width=9, height=6, units="in", res=175, quality=100)
ggplot(ld4_6, aes(log10(QVAL), fill=cytoband2)) + geom_histogram(binwidth = 0.05) + scale_fill_discrete(name="cytoband")
dev.off()


#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------- 170417 - CHROMOTHRIPSIS - PATHWAY AND GENE FREQUENCY/PERCENTAGE IN CASES AND CTRLS ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
library(reactome.db)

# take Rmatrix/LRmatrix (resTable) and a pathwayID:symbol table, find genes contained in each pathway and
# repeat each row as many times as genes are contained in each pathway
# column name for pathway id must be "pathID"
# column name for gene symbol must be "symbol"
getGeneSymbolsAndCreateNewRmatrix <- function(resTable, path2symbol){
  
  resTable2=NULL
  
  # loop over resTable rows (pathways)
  for(r in 1:dim(resTable)[1]){
    
    if(r%%10==0){cat("pathway", r, "/", dim(resTable)[1], "\n")}
    # get pathway[r]
    pathway = as.character(resTable$pathID[r])
    # get genes in pathway[r]
    genes = path2symbol[path2symbol$pathID==pathway,"symbol"]
    genes = genes[complete.cases(genes)]
    
    if(length(genes)>0){
      # repeat row [r] in resTable for N=length(genes) times and append everything to resTable2
      df=NULL
      for(l in 1:length(genes)){df = rbind(df, resTable[r,])}
      df$symbol = genes
      resTable2 = rbind(resTable2, df)
    }
    
  }
  
  return(resTable2)
  
}

# load LR results and fisher on genes results (loss/del)
logRegr = read.delim(paste0(amlDIR, "/chromothripsis/170417_logisticRegressionREACTOME_loss-del.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
lossDel = read.delim(paste0(amlDIR, "/chromothripsis/170414_geneFisher_lossDel.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
# assign gene names as rownames
rownames(lossDel) = lossDel$GENE

# create pathwayID:symbol table for REACTOME
reactome2symbol = select(reactome.db, keys=as.character(logRegr$pathID), column="ENTREZID", keytype = "PATHID") # output is a dataframe with entrez ids in the second column
reactome2symbol$SYMBOL = select(org.Hs.eg.db, keys=reactome2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
reactome2symbol = reactome2symbol[complete.cases(reactome2symbol),]
colnames(reactome2symbol) = c("pathID", "symbol")

# expand LRmatrix to have a row for each gene in each pathway for REACTOME
logRegr = getGeneSymbolsAndCreateNewRmatrix(logRegr, reactome2symbol)

# get values from geneFisher table into LRmatrix
logRegr$case = lossDel[logRegr$gene, "%Case"]
logRegr$ctrl = lossDel[logRegr$gene, "%Ctrl"]
logRegr$pval = lossDel[logRegr$gene, "PVAL"]
logRegr$qval = lossDel[logRegr$gene, "QVAL"]

# remove NAs
logRegr = logRegr[complete.cases(logRegr),]

# edit table
logRegr = logRegr[,c(1,2,3,6,7,8,9,10,11,12)]
colnames(logRegr)[c(3,4,5,7,8,9,10)] = c("Logistic Regression Coefficient", "pval", "Adjusted p-value", "%case", "%ctrl", "pval", "Adjusted p-value")

write.table(logRegr, file=paste0(amlDIR, "/chromothripsis/170417_path_gene_percentage_loss-del.tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)


#-------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------- 170414 - CHROMOTHRIPSIS - RANKING OF CHROMOSOMES FROM FISHER ON GENES -----------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
library(org.Hs.eg.db)
library(ggplot2)

# annotate gene coordinates for table containing gene symbols
# "gene_col" is the number of the column containing the gene symbols
getCoordinatesForSymbols <- function(genes, gene_col){
  
  # get gene list
  symbol = genes[,gene_col]
  
  # get symbol2entrezID table
  genes2entrez = select(org.Hs.eg.db, keys=as.character(symbol), columns="ENTREZID", keytype = "SYMBOL")
  colnames(genes2entrez)[2] = "entrez"
  
  # as it has multiple entrezID per symbol, get those symbols
  sel = genes2entrez[which(duplicated(genes2entrez$SYMBOL)), "SYMBOL"]
  dup = genes[genes[,gene_col] %in% sel,]
  # then add rows for duplicated symbols to "genes"
  genes = rbind(genes, dup)
  
  # sort both tables by symbol and merge
  genes = genes[order(genes[,gene_col]),]
  genes2entrez = genes2entrez[order(genes2entrez$SYMBOL),]
  genes = cbind(genes, entrez=genes2entrez$entrez, stringsAsFactors=F)
  
  # remove symbols with unknown entrezID
  genes = genes[complete.cases(genes),]
  
  # add columns
  genes$chr = rep(NA, dim(genes)[1])
  genes$start = rep(NA, dim(genes)[1])
  genes$end = rep(NA, dim(genes)[1])
  genes$strand = rep("+", dim(genes)[1])
  genes$chrArm = rep(NA, dim(genes)[1])
  
  # get list of start and end positions of entrez ids
  startPos = mget(genes$entrez,envir=org.Hs.egCHRLOC, ifnotfound=NA)
  endPos = mget(genes$entrez,envir=org.Hs.egCHRLOCEND, ifnotfound=NA)
  # get list of cytobands of entrez ids
  cytoband = mget(genes$entrez,envir=org.Hs.egMAP, ifnotfound=NA)
  
  # loop over table rows and fill start positions, end positions and strand sign of genes
  for(r in 1:dim(genes)[1]){
    if(r%%1000==0){cat("filling info for gene", r, "/", dim(genes)[1], "\n")}
    entrez = genes$entrez[r]
    if(!is.na(entrez)){
      if(!is.na(startPos[[entrez]][1])){
        genes$chr[r] = as.numeric(names(startPos[[entrez]][1]))
        genes$start[r] = as.character(startPos[[entrez]][1])
        genes$end[r] = as.character(endPos[[entrez]][1])
        genes$chrArm[r] = cytoband[[entrez]]
        if(grepl("-", genes$start[r])){
          genes$start[r] = as.numeric(strsplit(genes$start[r], "-")[[1]][2]) # positions on negative strands come as negative integers
          genes$end[r] = as.numeric(strsplit(genes$end[r], "-")[[1]][2]) # positions on negative strands come as negative integers
          genes$strand[r] = "-"}
      }
    }
  }
  
  #q = as.vector(sapply(genes[grep("q", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "q")[[1]][1], "q")))
  
  # extract chromosome arm info
  genes[grep("q", genes$chrArm), "chrArm"] = as.vector(sapply(genes[grep("q", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "q")[[1]][1], "q")))
  genes[grep("p", genes$chrArm), "chrArm"] = as.vector(sapply(genes[grep("p", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "p")[[1]][1], "p")))
  
  # remove symbols with unknown chromosome
  genes = genes[!is.na(genes$chrArm),]
  
  genes$start = as.numeric(genes$start)
  genes$end = as.numeric(genes$end)
  
  return(genes)
  
}

# load sample list and select case and ctrl
sampleID = read.delim(paste0(amlDIR, "/chromothripsis/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
selCase = sampleID[sampleID$chromothripsis==1,] # get case samples
selCtrl = sampleID[sampleID$chromothripsis==2,] # get ctrl samples

# load results of fisher on genes for loss/del
lossDel = read.delim(paste0(amlDIR, "/chromothripsis/170414_geneFisher_lossDel.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)

# add entrezID and chromosome coordinates
lossDel = getCoordinatesForSymbols(lossDel, 1)

# get chrArm names
chrArms = names(table(lossDel$chrArm))

# make table with summary() for qvals in each chromosomes
qvals = NULL
for(a in 1:length(chrArms)){
  summary = summary(lossDel[lossDel$chrArm == chrArms[a], "QVAL"])
  qvals = as.data.frame(rbind(qvals, c(chrArms[a], as.numeric(summary))), stringsAsFactors=F)}

# convert to numeric, assign colnames and sort by medians
qvals = data.frame(qvals[,1], apply(qvals[-1], 2, as.numeric), stringsAsFactors = F)
colnames(qvals) = c("chrArm", "min", "1st_qu", "median", "mean", "3rd_qu", "max")
qvals = qvals[order(qvals$median),]

#write.table(qvals, file=paste0(amlDIR, "/chromothripsis/170417_median_qval_chromosomes_loss-del.tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(qvals, file=paste0(amlDIR, "/chromothripsis/170417_median_qval_chromosomeArms_loss-del.tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- 170414 - CHROMOTHRIPSIS - FISHER ON CHROM vs TP53 AMONG 5qloss ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

library(IRanges)
library(GenomicRanges)

# create "segments" table with the coordinates of all altered segments and convert it to GRanges
# "table" must be an edited output from NEXUS, "order" must be a table specifing the radial position of each patient
getGRangesFromNexus <- function(table){
  
  gainDup = NULL
  lossDel = NULL
  
  for(r in 1:dim(table)[1]){
    
    # get event type
    event = table$Event[r]
    
    # get genomic coordinates
    coordinates = table$Chromosome.Region[r]
    
    # extract chr number and save as "hs#"
    chr = strsplit(coordinates, ":")[[1]][1]
    
    # extract "start" and "end" positions and convert to format "xxxxxxx"
    positions = strsplit(coordinates, ":")[[1]][2]
    
    start = strsplit(positions, "-")[[1]][1]
    start = as.numeric(gsub(",", "", start))
    
    end = strsplit(positions, "-")[[1]][2]
    end = as.numeric(gsub(",", "", end))
    
    patient = table$Sample[r]

    # fill info for GRanges of a specific event
    if(event == "gain" | event == "dup"){gainDup = rbind(gainDup, data.frame(chr=chr, start=start, end=end, patient=patient, stringsAsFactors = F))}
    if(event == "loss" | event == "del"){lossDel = rbind(lossDel, data.frame(chr=chr, start=start, end=end, patient=patient, stringsAsFactors = F))}
    
  }
  
  gainDup = makeGRangesFromDataFrame(gainDup, keep.extra.columns = T)
  lossDel = makeGRangesFromDataFrame(lossDel, keep.extra.columns = T)
  
  return(list(gainDup, lossDel))
  
}


# load and edit sample table
sampleID = read.delim(paste0(amlDIR, "/chromothripsis/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
selCase = sampleID[sampleID$chromothripsis==1,] # get case samples
selCtrl = sampleID[sampleID$chromothripsis==2,] # get ctrl samples
totCase = dim(selCase)[1] # get number of cases
totCtrl = dim(selCtrl)[1] # get number of ctrls

# load tables for segments
dgv_filtered = read.delim(paste0(amlDIR, "/chromothripsis/dgv_filtered.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
dgv_filtered$Sample<-unlist(lapply(strsplit(dgv_filtered$Sample, " "), "[[", 1))

# remove undesired samples
segments = dgv_filtered[which(dgv_filtered$Sample %in% sampleID$Sample),]

# remove LOH events and change other event names
segments = segments[segments$Event!="LOH",]
segments[segments$Event=="CN Loss", "Event"] = "loss"
segments[segments$Event=="CN Gain", "Event"] = "gain"
segments[segments$Event=="High Copy Gain", "Event"] = "dup"
segments[segments$Event=="Homozygous Copy Loss", "Event"] = "del"

segmentsCase = segments[segments$Sample %in% selCase$Sample,]
segmentsCtrl = segments[segments$Sample %in% selCtrl$Sample,]

# extract genomic coordinate data and build GRanges for gainDup and lossDel
case = getGRangesFromNexus(segmentsCase)
ctrl = getGRangesFromNexus(segmentsCtrl)

gdCase = case[[1]]
ldCase = case[[2]]
gdCtrl = ctrl[[1]]
ldCtrl = ctrl[[2]]

chr5case = ldCase[which(ldCase@seqnames == "chr5")]
chr5ctrl = ldCtrl[which(ldCtrl@seqnames == "chr5")]

cytoband = read.delim(paste0(circosTutorialDIR, "/karyotype.human.txt"), header=T, sep=" ", stringsAsFactors = F, check.names = F)
band5 = cytoband[cytoband$hs1="hs5",]

ds = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis/170414_5qloss_tp53.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "na")
ds$Sample = sapply(ds$Sample, function(x) strsplit(x," ")[[1]][1])

# FISHER_1. sample=ALL; case/ctrl=chromothripsis; pos/neg=5qloss
selCase = ds[ds$chromothripsis=="yes",]
selCtrl = ds[ds$chromothripsis=="no",]
totCase = dim(selCase)[1]
totCtrl = dim(selCtrl)[1]
posCase = dim(selCase[selCase$`5q31.1_5q31.3_loss`=="yes",])[1]
posCtrl = dim(selCtrl[selCtrl$`5q31.1_5q31.3_loss`=="yes",])[1]
negCase = totCase-posCase
negCtrl = totCtrl-posCtrl

x1 = matrix(c(posCase, negCase, posCtrl, negCtrl), nrow = 2, ncol=2)
f1 = fisher.test(x1, alternative = "two.sided", conf.int = T, conf.level = 0.99)
p1 = f1$p.value
est1 = f1$estimate

# FISHER_2. sample=5qloss; case/ctrl=chromothripsis; pos/neg=TP53mut
ds2 = ds[ds$`5q31.1_5q31.3_loss`=="yes",]

selCase = ds2[ds2$chromothripsis=="yes",]
selCtrl = ds2[ds2$chromothripsis=="no",]
totCase = dim(selCase)[1]
totCtrl = dim(selCtrl)[1]
posCase = dim(selCase[selCase$TP53_mut=="yes",])[1]
posCtrl = dim(selCtrl[selCtrl$TP53_mut=="yes",])[1]
negCase = totCase-posCase
negCtrl = totCtrl-posCtrl

x2 = matrix(c(posCase, negCase, posCtrl, negCtrl), nrow = 2, ncol=2)
f2 = fisher.test(x2, alternative = "two.sided", conf.int = T, conf.level = 0.99)
p2 = f2$p.value
est2 = f2$estimate

# FISHER_3. sample=NON_5qloss; case/ctrl=chromothripsis; pos/neg=TP53mut
ds3 = ds[ds$`5q31.1_5q31.3_loss`=="no",]

selCase = ds3[ds3$chromothripsis=="yes",]
selCtrl = ds3[ds3$chromothripsis=="no",]
totCase = dim(selCase)[1]
totCtrl = dim(selCtrl)[1]
posCase = dim(selCase[selCase$TP53_mut=="yes",])[1]
posCtrl = dim(selCtrl[selCtrl$TP53_mut=="yes",])[1]
negCase = totCase-posCase
negCtrl = totCtrl-posCtrl

x3 = matrix(c(posCase, negCase, posCtrl, negCtrl), nrow = 2, ncol=2)
f3 = fisher.test(x3, alternative = "two.sided", conf.int = T, conf.level = 0.99)
p3 = f3$p.value
est3 = f3$estimate






#-------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- 170404 - ANEUPLOIDY - CHECK IF MUTATED GENES ALSO HAVE CNAs ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------


gainDup = read.delim(paste0(amlDIR, "/aneuploidy/170403_gain-dup_toCompareToWES.txt"), header=F, sep="\t", stringsAsFactors = F, na.strings="")
colnames(gainDup) = gainDup[2,]
rownames(gainDup)[-1] = gainDup[-1,2]

lossDel = read.delim(paste0(amlDIR, "/aneuploidy/170403_loss-del_toCompareToWES.txt"), header=F, sep="\t", stringsAsFactors = F, na.strings="")
colnames(lossDel) = lossDel[1,]
rownames(lossDel) = lossDel[,2]

mut = read.delim(paste0(amlDIR, "/aneuploidy/170403_mut_toCompareToSNP.txt"), header=T, sep="\t", stringsAsFactors = F, na.strings="")

gdGenes = unique(colnames(gainDup)[-(1:2)])
gdPatients = unique(rownames(gainDup)[-(1:2)])

ldGenes = unique(colnames(lossDel)[-(1:2)])
ldPatients = unique(rownames(lossDel)[-1])

for(r in 1:dim(mut)[1]){
  print(r)
  gene = mut$GENE[r]
  patient = mut$Patient.ID[r]
  if(patient %in% gdPatients){
    if(gene %in% gdGenes){
      if(gainDup[patient,gene]>0){mut$gain.dup[r] = "y"
      }else{mut$gain.dup[r] = "n"}
    }else{mut$gain.dup[r] = "n"}
  
  if(patient %in% ldPatients){
    if(gene %in% ldGenes){
      if(lossDel[patient,gene]>0){mut$loss.del[r] = "y"
      }else{mut$loss.del[r] = "n"}
    }else{mut$loss.del[r] = "n"}}}}

write.table(mut, file=paste0(amlDIR, "/aneuploidy/170403_mut_comparedToSNP.txt"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)




#-------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- 170403 - ANEUPLOIDY - PATHWAY AND GENE FREQUENCY/PERCENTAGE IN CASES AND CTRLS ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# take a string "genes" like --> "gene1 (coefficient1); gene2 (coefficient2); .....; geneN (coefficientN)"
# extract all genes
# look up updated coefficients from "percentageTable"
# and re-concatenate all genes and coefficients like -->  "gene1 (NEWcoefficient1); gene2 (NEWcoefficient2); .....; geneN (NEWcoefficientN)"
undoAndRedoGeneWithPercentage <- function(genes, percentageTable){
  
  # split string by "; "
  genes = unlist(strsplit(genes, "; "))
  
  # get only gene names (i.e. remove " (coefficient)" from all resulting strings) and convert to data frame
  genes = data.frame(gene=unlist(lapply(strsplit(genes, " "), "[[", 1)), stringsAsFactors = F)
  
  # add column for NEW coefficients
  genes$percentage = rep(0, dim(genes)[1])
  
  for(r2 in 1:dim(genes)[1]){
    gene = genes$gene[r2]
    genes$percentage[r2] = percentageTable[percentageTable$gene == gene, 10][1]}
  
  genes$merge = paste0(genes[,1], " (", genes[,2], ");")
  
  genes = paste(genes$merge, collapse = " ")
  
  return(genes)
  
}


gainDup = read.delim(paste0(amlDIR, "/aneuploidy/170329_pathwayFisher_genes_GOBP_gain-dup.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
enrichTable = read.delim(paste0(amlDIR, "/aneuploidy/170403_enrichment_table_to_update.txt"), header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = c("","#N/A"))
logRegr = read.delim(paste0(amlDIR, "/aneuploidy/170403_logisticRegressionGOBP_gain-dup.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)

colnames(enrichTable) = c("pathID", "pathName", "path_case", "path_ctrl", "qval", "gene_case", "gene_ctrl")

# round to 1 digits values in columns 3,4,10,11 of gainDup
gainDup[,c(3,4,10,11)] = round(gainDup[,c(3,4,10,11)], digits=1)
# round to 4 digits values in "qval" column of enrichTable
enrichTable$qval = round(enrichTable$qval, digits=4)

# loop over pathways (rows)
for(r in 1:dim(enrichTable)[1]){
  
  # if it's not a line with a title
  if(!is.na(enrichTable$pathName[r])){
    
    # get pathway id
    pathway = enrichTable$pathID[r]
    
    # check in gainDup and update percentage of cases and ctrls with altered pathway
    enrichTable[r,c(3,4)] = gainDup[gainDup$PATHID == pathway, c(3,4)][1,] #[1,] is because in gainDup table each pathway appear in several rows
    
    # get string with all genes and gene percentages in pathway for cases
    gene_case = enrichTable$gene_case[r]
    # for cases, split string in single genes, update coefficients and re-collate
    if(!is.na(gene_case)){enrichTable$gene_case[r] = undoAndRedoGeneWithPercentage(gene_case, gainDup)}
    
    # get string with all genes and gene percentages in pathway for ctrls
    gene_ctrl = enrichTable$gene_ctrl[r]
    # for ctrls, split string in single genes, update coefficients and re-collate
    if(!is.na(gene_ctrl)){enrichTable$gene_ctrl[r] = undoAndRedoGeneWithPercentage(gene_ctrl, gainDup)}
    
  }
}

# assign pathID as rownames
rownames(logRegr) = logRegr$pathID

# round up logistic regression coefficients
enrichTable$Estimate = round(logRegr[enrichTable$pathID, "Estimate"], digits = 3)
#enrichTable$qval = round(enrichTable$qval, digits=4)

# select desired columns and edit names
enrichTable = enrichTable[,c(1,2,8,5,3,4,6,7)]
colnames(enrichTable)[3:4] = c("Logistic regression coefficient", "P-value (adj)")

write.table(enrichTable, file="/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/170403_enrichment_table_updated.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F, na="")



#-------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- 170329 - ANEUPLOIDY - PATHWAY AND GENE FREQUENCY/PERCENTAGE IN CASES AND CTRLS ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# take a string "genes" like --> "gene1 (coefficient1); gene2 (coefficient2); .....; geneN (coefficientN)"
# extract all genes
# look up updated coefficients from "percentageTable"
# and re-concatenate all genes and coefficients like -->  "gene1 (NEWcoefficient1); gene2 (NEWcoefficient2); .....; geneN (NEWcoefficientN)"
undoAndRedoGeneWithPercentage <- function(genes, percentageTable){
  
  # split string by "; "
  genes = unlist(strsplit(genes, "; "))
  
  # get only gene names (i.e. remove " (coefficient)" from all resulting strings) and convert to data frame
  genes = data.frame(gene=unlist(lapply(strsplit(genes, " "), "[[", 1)), stringsAsFactors = F)
  
  # add column for NEW coefficients
  genes$percentage = rep(0, dim(genes)[1])
  
  for(r2 in 1:dim(genes)[1]){
    gene = genes$gene[r2]
    genes$percentage[r2] = percentageTable[percentageTable$gene == gene, 10][1]}
  
  genes$merge = paste0(genes[,1], " (", genes[,2], ");")
  
  genes = paste(genes$merge, collapse = " ")
  
  return(genes)
  
}


gainDup = read.delim(paste0(amlDIR, "/aneuploidy/170329_pathwayFisher_genes_GOBP_gain-dup.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
enrichTable = read.delim(paste0(amlDIR, "/aneuploidy/170329_enrichment_table_to_update.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "")

colnames(enrichTable) = c("pathID", "pathName", "path_case", "path_ctrl", "qval", "gene_case", "gene_ctrl")

# round to 1 digits values in columns 3,4,10,11 of gainDup
gainDup[,c(3,4,10,11)] = round(gainDup[,c(3,4,10,11)], digits=1)
# round to 4 digits values in "qval" column of enrichTable
enrichTable$qval = round(enrichTable$qval, digits=4)

# loop over pathways (rows)
for(r in 1:dim(enrichTable)[1]){
  
  # if it's not a line with a title
  if(!is.na(enrichTable$pathName[r])){
    
    # get pathway id
    pathway = enrichTable$pathID[r]
    
    # check in gainDup and update percentage of cases and ctrls with altered pathway
    enrichTable[r,c(3,4)] = gainDup[gainDup$PATHID == pathway, c(3,4)][1,] #[1,] is because in gainDup table each pathway appear in several rows
    
    # get string with all genes and gene percentages in pathway for cases
    gene_case = enrichTable$gene_case[r]
    # for cases, split string in single genes, update coefficients and re-collate
    if(!is.na(gene_case)){enrichTable$gene_case[r] = undoAndRedoGeneWithPercentage(gene_case, gainDup)}
    
    # get string with all genes and gene percentages in pathway for ctrls
    gene_ctrl = enrichTable$gene_ctrl[r]
    # for ctrls, split string in single genes, update coefficients and re-collate
    if(!is.na(gene_ctrl)){enrichTable$gene_ctrl[r] = undoAndRedoGeneWithPercentage(gene_ctrl, gainDup)}
    
  }
}

write.table(enrichTable, file="/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/170329_enrichment_table_updated.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F, na="")





#----------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- JANUARY 2017 - CHROMOTHRIPSIS - DATASET PRE-PROCESSING ----------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis")

# # NOT REQUIRED ANYMORE
# # create list of genes on chrX and Y from a GENCODE annotation list from Table Browser
# gencode = read.table("gencode.txt", header = T, sep = "\t", quote="", stringsAsFactors = F)
# chrxy = gencode[(gencode[,3]=="chrX") | (gencode[,3]=="chrY") , c(3,13)]
# colnames(chrxy) = c("chr", "gene")
# write.table(chrxy, file="genes_chrXY.tsv", sep="\t", dec=".", row.names=F, quote=F)

# create unique list of all samples in two NEXUS output (kralovic and malek)
kralovic = read.delim("kralovic.tsv", header = T, sep = "\t", quote="", stringsAsFactors = F)
malek = read.delim("malek.tsv", header = T, sep = "\t", quote="", stringsAsFactors = F)
samples = sort(unique(c(unique(kralovic$Sample), unique(malek$Sample))))
write.table(samples, file="samples_list.tsv", sep="\t", col.names = F, row.names=F, dec=".", quote=F)

# adapt columns and merge "kralovic" and "malek" into one table
sel = which(colnames(kralovic) %in% colnames(malek))
kralovic = kralovic[, sel]
dgv_filtered = rbind(kralovic, malek)
write.table(dgv_filtered, file="dgv_filtered.tsv", sep="\t", col.names = T, row.names=F, dec=".", quote=F)

# compare samples list obtained from Marco's files to samples list from MCF
list1 = read.delim("samples_list.tsv", header = F, sep = "\t", quote="", stringsAsFactors = F)
list2 = read.delim("case_crtl_170120.tsv", header = T, sep = "\t", quote="", stringsAsFactors = F)
sel = which(list2[,1] %in% list1[,1])
unmatched = list2[-sel,]



#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------- DATE UNKNOWN - ANEUPLOIDY - ADD GENE SYMBOL TO REACTOME ANNOTATION FILE ----------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

library(org.Hs.eg.db)
library(annotate)

setwd("/media/DATA1/aml/ProjectAE/last/tables_revised/")

reactome = list.files(getwd(), "AnnoReactome")

one = read.table(reactome[1], header=T, sep="\t", quote="", stringsAsFactors=F)
two = data.frame(one[,1], vector(mode="character", length=dim(one)[1]), one[,2], one[,3])
colnames(two) = c(names(one)[1], "gene", names(one)[2], names(one)[3]) 

for(i in 1:dim(two)[1]){
  two$gene = lookUp(as.character(two$ENTREZID), "org.Hs.eg.db", "SYMBOL")
}






#-------------------- test ----------------


