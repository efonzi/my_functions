
library(org.Hs.eg.db)
library(plyr)

library(IRanges)
library(GenomicRanges)

testDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/circos/test"
aneuplDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/circos/aneupl"
chromotDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/circos/chromot"
profDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/circos/prof"

#----------------------------------------------------------------- FUNCTIONS ----------------------------------------------------------------------
# annotate gene coordinates from gene symbol list
getCoordinatesForSymbols <- function(genes){
  
  # transform to entrez ids
  genes$entrez = select(org.Hs.eg.db, keys=as.character(genes$symbol), columns="ENTREZID", keytype = "SYMBOL")[,2]
  
  # add empty columns
  genes$chr = rep(0, dim(genes)[1])
  genes$start = rep(0, dim(genes)[1])
  genes$end = rep(0, dim(genes)[1])
  genes$strand = rep("+", dim(genes)[1])
  
  # get list of start and end positions of entrez ids
  startPos = mget(genes$entrez,envir=org.Hs.egCHRLOC, ifnotfound=NA)
  endPos = mget(genes$entrez,envir=org.Hs.egCHRLOCEND, ifnotfound=NA)
  
  # loop over table rows and fill start positions, end positions and strand sign of genes
  for(r in 1:dim(genes)[1]){
    genes$chr[r] = as.numeric(names(startPos[[r]][1]))
    genes$start[r] = as.character(startPos[[r]][1])
    genes$end[r] = as.character(endPos[[r]][1])
    if(grepl("-", genes$start[r])){
      genes$start[r] = as.numeric(strsplit(genes$start[r], "-")[[1]][2]) # positions on negative strands come as negative integers
      genes$end[r] = as.numeric(strsplit(genes$end[r], "-")[[1]][2]) # positions on negative strands come as negative integers
      genes$strand[r] = "-"}
  }
  
  genes$start = as.numeric(genes$start)
  genes$end = as.numeric(genes$end)
  
  return(genes)
  
}

# create "segments" table with the coordinates of all altered segments
# "table" must be an edited output from NEXUS, "order" must be a table specifing the radial position of each patient and
# "palette" must be a table assigning color to specific events
createSegmentsForHighlight <- function(table, order, palette){
  
  l = dim(table)[1]
  
  highlight = data.frame(chr=rep(".",l), start=rep(0, l), end=rep(0, l), options=rep(".", l), stringsAsFactors = F)
  
  gainDup = NULL
  lossDel = NULL
  
  for(r in 1:l){
    
    # get event type
    event = table$Event[r]
    
    # get coordinates from table in the format "chr#:x,xxx,xxx-x,xxx,xxx"
    coordinates = table$Chromosome.Region[r]
    
    # extract chr number and save as "hs#"
    chr = strsplit(coordinates, ":")[[1]][1]
    chr = paste0("hs", strsplit(chr, "chr")[[1]][2])
    
    # extract "start" and "end" positions and convert to format "xxxxxxx"
    positions = strsplit(coordinates, ":")[[1]][2]
    
    start = strsplit(positions, "-")[[1]][1]
    start = as.numeric(gsub(",", "", start))
    
    end = strsplit(positions, "-")[[1]][2]
    end = as.numeric(gsub(",", "", end))
    
    # determine color of event
    color = palette[palette$event==event, "color"]
    
    # get radius coordinates (r0 and r1) for patient
    patient = table$Sample[r]
    r0 = order[order$Sample==patient, "r0"]
    r1 = order[order$Sample==patient, "r1"]
    
    # fill info for "highlight" dataset --> chromosome, start, end and options(color, r0, r1)
    highlight$chr[r] = chr
    highlight$start[r] = start
    highlight$end[r] = end
    highlight$options[r] = paste0("fill_color=", color, ",r0=", r0, "r,r1=", r1, "r")
    
    # fill info for GRanges of a specific event
    if(event == "gain" | event == "dup"){gainDup = rbind(gainDup, data.frame(chr=chr, start=start, end=end, stringsAsFactors = F))}
    if(event == "loss" | event == "del"){lossDel = rbind(lossDel, data.frame(chr=chr, start=start, end=end, stringsAsFactors = F))}

  }
  
  gainDup = makeGRangesFromDataFrame(gainDup)
  lossDel = makeGRangesFromDataFrame(lossDel)
  
  return(list(highlight, gainDup, lossDel))
  
}

# create co-occurrence table for LINKS --> a table with 1:1 gene combinations with counts of patients that harbor them
# PATIENTxGENE matrix for gainDup and for lossDel are necessary (any event could be suitable)
# column with gene names in "genes" must be labelled "symbol"
computeCombinationsForLinks <- function(gainDup, lossDel, genes){
  
  # get only desired genes (columns)
  links_gd = gainDup[,which(colnames(gainDup) %in% genes$symbol)]
  links_ld = lossDel[,which(colnames(lossDel) %in% genes$symbol)]
  
  # merge gain-dup table to loss-del table
  links_gdld = cbind(links_gd, links_ld)
  
  # substitute events with corresponding gene name, in order to compute combinations
  for(c in 1:dim(links_gdld)[2]){
    gene = colnames(links_gdld)[c]
    links_gdld[links_gdld[,c]>0,c] = gene}
  
  # for each patient calculate all 1:1 combinations of any type of event and create a total list
  combinations=NULL
  for(r in 1:dim(links_gdld)[1]){
    events = as.character(links_gdld[r,which(links_gdld[r,]!=0)])
    if(length(events)>1){
      df = as.data.frame(t(combn(events, 2)), stringsAsFactors=F) # combn() calculates all possible combinations of 2 events; t() transposes the resulting matrix
      combinations = rbind(combinations, df)}}
  
  colnames(combinations) = c("a","b")
  
  # count how many times each combination is repeated
  combinations = ddply(combinations,.(a,b),nrow)
 
  return(combinations) 
}

# create dataframes for links in A and links in E, with following structure:
# chr1 - start1 - end1 - chr2 - start2 - end2 - color=....
createLinksTable <- function(combinations, genes, paletteName, paletteTot, paletteType, max){
  
  # create paletteString
  paletteString = paste0("color=", paletteName, "-", paletteTot, "-", paletteType, "-")
  
  # get total number of combinations
  l = dim(combinations)[1]
  
  # create empty dataframe with necessary fields
  links = data.frame(chr1=rep(".",l), start1=rep(0,l), end1=rep(0,l), chr2=rep(".",l), start2=rep(0,l), end2=rep(0,l), options=rep(0,l), stringsAsFactors = F)
  
  # loop over rows and fill each cell
  for(r in 1:l){
    gene1 = combinations$a[r]
    gene2 = combinations$b[r]
    count = combinations$V1[r]
    links$chr1[r] = genes[genes$symbol==gene1, "chr"]
    links$start1[r] = genes[genes$symbol==gene1, "start"]
    links$end1[r] = genes[genes$symbol==gene1, "end2"]
    links$chr2[r] = genes[genes$symbol==gene2, "chr"]
    links$start2[r] = genes[genes$symbol==gene2, "start"]
    links$end2[r] = genes[genes$symbol==gene2, "end2"]
    
    # # assign color according to level of link "connection"
    # if(count < max){links$options[r] = paste0(paletteString, ((paletteTot/max)*count)%/%1+1)
    # }else{links$color[r] = paste0(paletteString, paletteTot)}}
    
    # to group the less "connected" links and assign them the same color (which is not the clearest one and it will be visible)
    if(count < max){num = ((paletteTot/max)*count)%/%1+1
      if(num > 2){links$options[r] = paste0(paletteString, num, ",z=", num)
      }else{links$options[r] = paste0(paletteString, "2,z=2")}
    }else{links$options[r] = paste0(paletteString, paletteTot, ",z=", paletteTot)}}
  
    # # to convert clear color coefficients to high color and viceversa (suitable only for scales from 1 to 11)
    # if(count < max){num = ((paletteTot/max)*count)%/%1+1
    #   links$options[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))
    # }else{num = paletteTot
    #   links$options[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))}}
  
  return(links)
  
}

#----------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------- EDIT HUMAN.KARYOTYPE.TXT ----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

karyo1 = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/circos/test/karyotype.human.txt", header=F, sep=" ", stringsAsFactors = F)

chr = karyo1[karyo1$V1 == "chr",]
chrGain = chr
chrLoss = chr
chrGain$V3 = paste0(chrGain$V3, "G")
chrLoss$V3 = paste0(chrLoss$V3, "L")

band = karyo1[karyo1$V1 == "band",]
bandGain = band
bandLoss = band
bandGain$V2 = paste0(bandGain$V2, "G")
bandLoss$V2 = paste0(bandLoss$V2, "L")

karyoGL = rbind(chrGain, chrLoss, bandGain, bandLoss)

write.table(karyoGL, file=paste0(aneuplDIR, "/karyoGainLoss.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)


#----------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------ ANEUPLOIDY ----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------- FOR WHOLE CHROMOSOMES - GAIN/LOSS SEPARATED - SINGLE HISTOGRAM --------------------------------------

# load PATIENTxGENE matrices and sample list
gainDup = read.delim(paste0(aneuplDIR, "/matrixPatientGene_aneupl_array_AE_gain-dup_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
lossDel = read.delim(paste0(aneuplDIR, "/matrixPatientGene_aneupl_array_AE_loss-del_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)

sampleID = read.delim(paste0(aneuplDIR, "/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
selCase = sampleID[sampleID$AE==1,] # get case samples
selCtrl = sampleID[sampleID$AE==2,] # get ctrl samples
totCase = dim(selCase)[1] # get number of cases
totCtrl = dim(selCtrl)[1] # get number of ctrls

# change "MLLT4" to "AFDN" gene name (it's the same gene, but there is no ENTREZ id for "MLLT4" version)
colnames(gainDup)[colnames(gainDup) == "MLLT4"] = "AFDN"
colnames(lossDel)[colnames(lossDel) == "MLLT4"] = "AFDN"

# change all values to 1
for(r in 1:dim(gainDup)[1]){
  gainDup[r,which(gainDup[r,]>1)[-1]] = 1
  lossDel[r,which(lossDel[r,]>1)[-1]] = 1}

# separate samples into cases and ctrls for PATIENTxGENE matrices
gdCase = gainDup[gainDup$Sample %in% selCase$Sample,]
gdCtrl = gainDup[gainDup$Sample %in% selCtrl$Sample,]
ldCase = lossDel[lossDel$Sample %in% selCase$Sample,]
ldCtrl = lossDel[lossDel$Sample %in% selCtrl$Sample,]

# load and edit genes list for aneuploidy
genes = read.delim(paste0(aneuplDIR, "/geneList.txt"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
genes = data.frame(symbol=genes$GENE, event=genes$`event type`, stringsAsFactors = F)
genes[genes$symbol == "MLLT4", "symbol"] = "AFDN"
genes[genes$symbol == "NPM1", "event"] = "loss/del"

# annotate coordinates
genes = getCoordinatesForSymbols(genes)
# adjust end position of all genes to make them visible on plot (start + 1000000)
genes$end2 = genes$start + 500000
# edit chr names to CIRCOS format
genes$chr = paste0("hs", genes$chr)
# add "L" or "G" to chr names (because each chr will have to show only loss/del or gain/dup events)
genes[genes$event == "loss/del","chr"] = sapply(genes[genes$event == "loss/del","chr"], function(x) x=paste0(x, "L"))
genes[genes$event == "gain/dup","chr"] = sapply(genes[genes$event == "gain/dup","chr"], function(x) x=paste0(x, "G"))

# count events for cases and ctrl and calculate percentage
genes$case = rep(0, dim(genes)[1])
genes$ctrl = rep(0, dim(genes)[1])
for(r in 1:dim(genes)[1]){
  gene = genes$symbol[r]
  event = genes$event[r]
  if(event == "loss/del"){if(gene %in% colnames(ldCase)){genes$case[r] = sum(ldCase[,gene])/totCase*100}}
  if(event == "loss/del"){if(gene %in% colnames(ldCtrl)){genes$ctrl[r] = sum(ldCtrl[,gene])/totCtrl*100}}
  if(event == "gain/dup"){if(gene %in% colnames(gdCase)){genes$case[r] = sum(gdCase[,gene])/totCase*100}}
  if(event == "gain/dup"){if(gene %in% colnames(gdCtrl)){genes$ctrl[r] = sum(gdCtrl[,gene])/totCtrl*100}}}

# add options columns with color specification for CIRCOS histograms and labels
genes$options = rep(".", dim(genes)[1])
genes[genes$event == "loss/del", "options"] = "vdgreen"
genes[genes$event == "gain/dup", "options"] = "vdred"

# create table with gene positions for histogram track ("hist.txt")
histCase = cbind(genes[,c(4,5,8,9)], options=paste0("fill_color=", genes$options), stringsAsFactors=F)
histCtrl = cbind(genes[,c(4,5,8,10)], options=paste0("fill_color=", genes$options), stringsAsFactors=F)

# create table with gene LABELS for text track ("text.txt")
text = cbind(genes[,c(4,5,8,1)], options=paste0("color=", genes$options), stringsAsFactors=F)


combinationsCase = computeCombinationsForLinks(gdCase, ldCase, genes)
combinationsCtrl = computeCombinationsForLinks(gdCtrl, ldCtrl, genes)

# get highest number of patients with the same combination, among both case and ctrl combinations
max = max(max(combinationsCase$V1), max(combinationsCtrl$V1))

paletteName = "blues"
#paletteName = "ylgnbu"
#paletteName = "rdylgn"

paletteTot = 9
#paletteTot = 11

paletteType = "seq"
#paletteType = "div"

# create dataframe formatted for LINKS
linksCase = createLinksTable(combinationsCase, genes, paletteName, paletteTot, paletteType, max)
linksCtrl = createLinksTable(combinationsCtrl, genes, paletteName, paletteTot, paletteType, max)


write.table(histCase, file=paste0(aneuplDIR, "/histCase.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histCtrl, file=paste0(aneuplDIR, "/histCtrl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(text, file=paste0(aneuplDIR, "/textSeparate.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)

write.table(linksCase, file=paste0(aneuplDIR, "/linksCaseSeparate.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(linksCtrl, file=paste0(aneuplDIR, "/linksCtrlSeparate.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)


#---------------------------------------- FOR WHOLE CHROMOSOMES - GAIN/LOSS TOGETHER - DOUBLE HISTOGRAM -------------------------------------

# load PATIENTxGENE matrices and sample list
gainDup = read.delim(paste0(aneuplDIR, "/matrixPatientGene_aneupl_array_AE_gain-dup_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
lossDel = read.delim(paste0(aneuplDIR, "/matrixPatientGene_aneupl_array_AE_loss-del_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)

sampleID = read.delim(paste0(aneuplDIR, "/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
selCase = sampleID[sampleID$AE==1,] # get case samples
selCtrl = sampleID[sampleID$AE==2,] # get ctrl samples
totCase = dim(selCase)[1] # get number of cases
totCtrl = dim(selCtrl)[1] # get number of ctrls

# change "MLLT4" to "AFDN" gene name (it's the same gene, but there is no ENTREZ id for "MLLT4" version)
colnames(gainDup)[colnames(gainDup) == "MLLT4"] = "AFDN"
colnames(lossDel)[colnames(lossDel) == "MLLT4"] = "AFDN"

# change all values to 1
for(r in 1:dim(gainDup)[1]){
  gainDup[r,which(gainDup[r,]>1)[-1]] = 1
  lossDel[r,which(lossDel[r,]>1)[-1]] = 1}

# separate samples into cases and ctrls for PATIENTxGENE matrices
gdCase = gainDup[gainDup$Sample %in% selCase$Sample,]
gdCtrl = gainDup[gainDup$Sample %in% selCtrl$Sample,]
ldCase = lossDel[lossDel$Sample %in% selCase$Sample,]
ldCtrl = lossDel[lossDel$Sample %in% selCtrl$Sample,]

# load and edit genes list for aneuploidy
genes = read.delim(paste0(aneuplDIR, "/geneList.txt"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
genes = data.frame(symbol=genes$GENE, event=genes$`event type`, stringsAsFactors = F)
genes[genes$symbol == "MLLT4", "symbol"] = "AFDN"
genes[genes$symbol == "NPM1", "event"] = "loss/del"

# annotate coordinates
genes = getCoordinatesForSymbols(genes)
# adjust end position of all genes to make them visible on plot (start + 1000000)
genes$end2 = genes$start + 500000
# edit chr names to CIRCOS format
genes$chr = paste0("hs", genes$chr)

# count events of gain+dup and loss+del for cases and ctrl
genes$gdCase = rep(0, dim(genes)[1])
genes$ldCase = rep(0, dim(genes)[1])
genes$gdCtrl = rep(0, dim(genes)[1])
genes$ldCtrl = rep(0, dim(genes)[1])
for(r in 1:dim(genes)[1]){
  gene = genes$symbol[r]
  if(gene %in% colnames(gdCase)){genes$gdCase[r] = sum(gdCase[,gene])/totCase*100}
  if(gene %in% colnames(ldCase)){genes$ldCase[r] = sum(ldCase[,gene])/totCase*100}
  if(gene %in% colnames(gdCtrl)){genes$gdCtrl[r] = sum(gdCtrl[,gene])/totCtrl*100}
  if(gene %in% colnames(ldCtrl)){genes$ldCtrl[r] = sum(ldCtrl[,gene])/totCtrl*100}}

# create table with gene positions for histogram track ("hist.txt")
histGDcase = genes[,c(4,5,8,9)]
histLDcase = genes[,c(4,5,8,10)]
histGDctrl = genes[,c(4,5,8,11)]
histLDctrl = genes[,c(4,5,8,12)]

# create table with gene LABELS for text track ("text.txt")
textALL = cbind(genes[,c(4,5,8,1)], options=rep(".", dim(genes)[1]), stringsAsFactors = F)
for(r in 1:dim(textALL)[1]){
  if(genes$event[r] == "loss/del"){textALL$options[r] = "color=vdgreen"}
  if(genes$event[r] == "gain/dup"){textALL$options[r] = "color=vdred"}}

textALL = cbind(genes[,c(4,5,8,1)], stringsAsFactors = F) # for prof circos

combinationsCase = computeCombinationsForLinks(gdCase, ldCase, genes)
combinationsCtrl = computeCombinationsForLinks(gdCtrl, ldCtrl, genes)

# get highest number of patients with the same combination, among both case and ctrl combinations
max = max(max(combinationsCase$V1), max(combinationsCtrl$V1))

#paletteName = "blues"
paletteName = "ylgnbu"
#paletteName = "rdylgn"

paletteTot = 9
#paletteTot = 11

paletteType = "seq"
#paletteType = "div"

# create dataframe formatted for LINKS
linksCase = createLinksTable(combinationsCase, genes, paletteName, paletteTot, paletteType, max)
linksCtrl = createLinksTable(combinationsCtrl, genes, paletteName, paletteTot, paletteType, max)

write.table(histGDcase, file=paste0(aneuplDIR, "/histGainDupAneupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histLDcase, file=paste0(aneuplDIR, "/histLossDelAneupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histGDctrl, file=paste0(aneuplDIR, "/histGainDupEupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histLDctrl, file=paste0(aneuplDIR, "/histLossDelEupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(textALL, file=paste0(aneuplDIR, "/textTogether.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)

write.table(linksCase, file=paste0(aneuplDIR, "/linksCaseTogether.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(linksCtrl, file=paste0(aneuplDIR, "/linksCtrlTogether.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)

write.table(linksCase, file=paste0(profDIR, "/linksCaseProf.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(textALL, file=paste0(profDIR, "/textProf.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)


#---------------------------------------- FOR CYTOBANDS --------------------------------------

# load human karyotype file
karyo = read.delim(paste0(testDIR, "/karyotype.human.txt"), header=F, sep=" ", stringsAsFactors = F)
# subset for cytoband info
band = karyo[-(1:24),-1]
band[,1] = substr(band[,1], 3, nchar(band[,1]))
band$cytoband = paste0(band$V2, band$V3)

# tetraploid patients
tetra = c("20150331_MG_AML_Cyto58_(CytoScanHD_Array)", "20150520_MG_AML_Cyto85_(CytoScanHD_Array)")

# load genes list for aneuploidy
table1 = read.delim(paste0(aneuplDIR, "/geneList.txt"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
# get gene symbols
table2 = data.frame(symbol=table1$GENE, stringsAsFactors = F)
# transform to entrez ids
table2$entrez = select(org.Hs.eg.db, keys=as.character(table2$symbol), columns="ENTREZID", keytype = "SYMBOL")[,2]

# add empty columns
table2$chr = rep(0, dim(table2)[1])
table2$start = rep(0, dim(table2)[1])
table2$end = rep(0, dim(table2)[1])
table2$strand = rep("+", dim(table2)[1])
table2$bandSTART = rep(".", dim(table2)[1])
table2$bandEND = rep(".", dim(table2)[1])
table2$double_band = rep(".", dim(table2)[1])
table2$start_rel = rep(0, dim(table2)[1])
table2$end_rel = rep(0, dim(table2)[1])

# get list of start and end positions of entrez ids
startPos = mget(table2$entrez,envir=org.Hs.egCHRLOC, ifnotfound=NA)
endPos = mget(table2$entrez,envir=org.Hs.egCHRLOCEND, ifnotfound=NA)

# loop over table rows and fill start positions, end positions and strand sign of genes
for(r in 1:dim(table2)[1]){
  table2$chr[r] = as.numeric(names(startPos[[r]][1]))
  table2$start[r] = as.character(startPos[[r]][1])
  table2$end[r] = as.character(endPos[[r]][1])
  if(grepl("-", table2$start[r])){
    table2$start[r] = strsplit(table2$start[r], "-")[[1]][2] # positions on negative strands come as negative integers
    table2$end[r] = strsplit(table2$end[r], "-")[[1]][2] # positions on negative strands come as negative integers
    table2$strand[r] = "-"}
}

table2$start = as.numeric(table2$start)
table2$end = as.numeric(table2$end)

# loop over table rows, subset "band" table, loop over its rows and get cytoband info
for(r1 in 1:dim(table2)[1]){
  # get chromosome and start position info
  chrom = table2$chr[r1]
  start = table2$start[r1]
  end = table2$end[r1]
  # subset "band" for entries belonging to that chromosome
  bandTemp = band[band[,1]==chrom,]
  # loop over "bandTemp" rows, look for correct cytoband, extract info and save it to table
  for(r2 in 1:dim(bandTemp)[1]){
    if((bandTemp[r2,4]<start) & (bandTemp[r2,5]>start)){table2$bandSTART[r1] = paste0(chrom, bandTemp[r2,2])}
    if((bandTemp[r2,4]<end) & (bandTemp[r2,5]>end)){table2$bandEND[r1] = paste0(chrom, bandTemp[r2,2])}
  }
  # if gene lies on one cytoband write "==", if it spans two cytobands write "!="
  if(table2$bandSTART[r1]==table2$bandEND[r1]){table2$double_band[r1]="=="
  }else{table2$double_band[r1]="!="}
}

to_split=NULL
to_split = table2[table2$double_band=="!=",]
table2[table2$double_band=="!=","start"] = 0
to_split$end = 0
table2 = rbind(table2, to_split)

for(r in 1:dim(table2)[1]){
  # select genes spanning two cytobands
  if(table2$double_band[r]=="!="){
    # get start position of 2nd band
    start_2nd_band = band[band$cytoband==table2$bandEND[r],4]
    if(table2$start[r]==0){
      # the start position of 2nd band is used as the position to split the gene
      table2$start[r]=start_2nd_band
      # adjust cytobands names
      table2$bandSTART[r] = table2$bandEND[r]}
    if(table2$end[r]==0){
      # the start position of 2nd band is used as the position to split the gene
      table2$end[r]=start_2nd_band
      # adjust cytobands names
      table2$bandEND[r]=table2$bandSTART[r]}
  }
}

# get list of unique bands involved
bandList = sort(unique(c(table2$bandSTART, table2$bandEND)))

# create data frame with all involved bands and info according to CIRCOS specification for karyotype files
l = length(bandList)
band_circos = data.frame(chr=rep("chr", l))
band_circos$hyphen = rep("-", l)
band_circos$ID = bandList
band_circos$LABEL = bandList
band_circos$START = rep(0, l)
band_circos$END = rep(0, l)
band_circos$COLOR = rep(".", l)
for(r in 1:dim(band_circos)[1]){
  cytoband = band_circos$ID[r]
  band_circos$START[r] = band[band$cytoband==cytoband,4]
  band_circos$END[r] = band[band$cytoband==cytoband,5]
  band_circos$COLOR[r] = band[band$cytoband==cytoband,6]}

# calculate relative START and END positions for genes (subtracting cytoband start position)
for(r in 1:dim(table2)[1]){
  # get start position of 1st cytoband (in most cases this will suffice)
  x = band_circos[band_circos$ID==table2$bandSTART[r],5]
  # get start position of 2nd cytoband
  x1 = band_circos[band_circos$ID==table2$bandEND[r],5]
  # subtracts relevant cytoband start position from gene coordinates
  if(table2$double_band[r] == "=="){
    table2$start_rel[r] = table2$start[r] - x +1 #+1 is to prevent some positions to become multiple of 10^6 and being written as xe+06
    table2$end_rel[r] = table2$end[r] - x +1
  }else{# in case of gene spanning two cytobands
    if(table2$start[r]!=x1){
      table2$start_rel[r] = table2$start[r] - x +1
      table2$end_rel[r] = table2$end[r] - x +1
    }else{
      table2$start_rel[r] = table2$start[r] - x1 +1
      table2$end_rel[r] = table2$end[r] - x1 +1
    }
  }
}

# calculating relative START and END positions for "band_circos"
band_circos$END = band_circos$END - band_circos$START +1 #+1 is to prevent some positions to become multiple of 10^6 and being written as xe+06

band_circos$START = 0

# split between "p" and "q" cytobands, add "chromosome" column and reorder by chromosome
p = band_circos[which(grepl("p", band_circos$ID)),]
q = band_circos[which(grepl("q", band_circos$ID)),]
p[, "chromosome"] = as.numeric(unlist(lapply(strsplit(p$ID, "p"), "[[", 1)))
q[, "chromosome"] = as.numeric(unlist(lapply(strsplit(q$ID, "q"), "[[", 1)))
band_circos = rbind(p, q)

band_circos = band_circos[order(band_circos$chromosome),]

# remove chromosome column
band_circos = band_circos[,-(dim(band_circos)[2])]

# load PATIENTxGENE matrices and sample list
gd = read.delim(paste0(aneuplDIR, "/matrixPatientGene_aneupl_array_AE_gain-dup_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
ld = read.delim(paste0(aneuplDIR, "/matrixPatientGene_aneupl_array_AE_loss-del_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID = read.delim(paste0(aneuplDIR, "/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))

# change all values to 1
for(r in 1:dim(gd)[1]){
  gd[r,which(gd[r,]>1)[-1]] = 1
  ld[r,which(ld[r,]>1)[-1]] = 1}

# # remove tetraploids
# gd = gd[!gd$Sample %in% tetra,]
# ld = ld[!ld$Sample %in% tetra,]

# separate samples into cases and ctrls
selCase = sampleID[sampleID$AE==1,]
selCtrl = sampleID[sampleID$AE==2,]
gdCase = gd[gd$Sample %in% selCase$Sample,]
gdCtrl = gd[gd$Sample %in% selCtrl$Sample,]
ldCase = ld[ld$Sample %in% selCase$Sample,]
ldCtrl = ld[ld$Sample %in% selCtrl$Sample,]

# get total number of cases and ctrls
totCase = dim(selCase)[1]
totCtrl = dim(selCtrl)[1]

# count events of gain+dup and loss+del for aneuploidy and euploidy and save to table2
table2$gdCase = rep(0, dim(table2)[1])
table2$ldCase = rep(0, dim(table2)[1])
table2$gdCtrl = rep(0, dim(table2)[1])
table2$ldCtrl = rep(0, dim(table2)[1])
for(r in 1:dim(table2)[1]){
  gene = table2$symbol[r]
  if(gene %in% colnames(gdCase)){table2$gdCase[r] = sum(gdCase[,gene])/totCase*100}
  if(gene %in% colnames(ldCase)){table2$ldCase[r] = sum(ldCase[,gene])/totCase*100}
  if(gene %in% colnames(gdCtrl)){table2$gdCtrl[r] = sum(gdCtrl[,gene])/totCtrl*100}
  if(gene %in% colnames(ldCtrl)){table2$ldCtrl[r] = sum(ldCtrl[,gene])/totCtrl*100}}

# create table with gene positions for histogram ("hist.txt")
histGDcase = table2[,c(7,10,11,12)]
histLDcase = table2[,c(7,10,11,13)]
histGDctrl = table2[,c(7,10,11,14)]
histLDctrl = table2[,c(7,10,11,15)]

#### create co-occurrence tables for LINKS
# get only desired genes (columns)
links_gd_A = gdCase[,which(colnames(gdCase) %in% table2$symbol)]
links_ld_A = ldCase[,which(colnames(ldCase) %in% table2$symbol)]

links_gd_E = gdCtrl[,which(colnames(gdCtrl) %in% table2$symbol)]
links_ld_E = ldCtrl[,which(colnames(ldCtrl) %in% table2$symbol)]

# merge gain-dup table to loss-del table
links_gdld_A = cbind(links_gd_A, links_ld_A)
links_gdld_E = cbind(links_gd_E, links_ld_E)

# substitute events with corresponding gene name, in order to compute combinations
for(c in 1:dim(links_gdld_A)[2]){
  gene = colnames(links_gdld_A)[c]
  links_gdld_A[links_gdld_A[,c]>0,c] = gene}

for(c in 1:dim(links_gdld_E)[2]){
  gene = colnames(links_gdld_E)[c]
  links_gdld_E[links_gdld_E[,c]>0,c] = gene}

# for each patient calculate all 1:1 combinations of any type of event and create a total list
combinations_A=NULL
combinations_E=NULL
for(r in 1:dim(links_gdld_A)[1]){
  events = as.character(links_gdld_A[r,which(links_gdld_A[r,]!=0)])
  if(length(events)>1){
    df = as.data.frame(t(combn(events, 2)), stringsAsFactors=F) # combn() calculates all possible combinations of 2 events; t() transposes the resulting matrix
    combinations_A = rbind(combinations_A, df)}}

for(r in 1:dim(links_gdld_E)[1]){
  events = as.character(links_gdld_E[r,which(links_gdld_E[r,]!=0)])
  if(length(events)>1){
    df = as.data.frame(t(combn(events, 2)), stringsAsFactors=F)
    combinations_E = rbind(combinations_E, df)}}

colnames(combinations_A) = c("a","b")
colnames(combinations_E) = c("a","b")

# count how many times each combination is repeated
combinations_A = ddply(combinations_A,.(a,b),nrow)
combinations_E = ddply(combinations_E,.(a,b),nrow)

# get highest number of patients with the same combination, among both case and ctrl combinations
max = max(max(combinations_A$V1), max(combinations_E$V1))

#paletteName = "blues"
paletteName = "ylgnbu"
#paletteName = "rdylgn"

paletteTot = 9
#paletteTot = 11

paletteType = "seq"
#paletteType = "div"

paletteString = paste0("color=", paletteName, "-", paletteTot, "-", paletteType, "-")

# create dataframes for links in A and links in E, with following structure:
# cytoband1 - start1 - end1 - cytoband2 - start2 - end2 - color=....
l = dim(combinations_A)[1]
links_A = data.frame(cytoband1=rep(".",l), start1=rep(0,l), end1=rep(0,l), cytoband2=rep(".",l), start2=rep(0,l), end2=rep(0,l), color=rep(0,l), stringsAsFactors = F)
for(r in 1:l){
  gene1 = combinations_A$a[r]
  gene2 = combinations_A$b[r]
  count = combinations_A$V1[r]
  links_A$cytoband1[r] = table2[table2$symbol==gene1, "bandSTART"]
  links_A$start1[r] = table2[table2$symbol==gene1, "start_rel"]
  links_A$end1[r] = table2[table2$symbol==gene1, "end_rel"]
  links_A$cytoband2[r] = table2[table2$symbol==gene2, "bandSTART"]
  links_A$start2[r] = table2[table2$symbol==gene2, "start_rel"]
  links_A$end2[r] = table2[table2$symbol==gene2, "end_rel"]
  
  # # assign color according to level of link "connection"
  # if(count < max){links_A$color[r] = paste0(paletteString, ((paletteTot/max)*count)%/%1+1)
  # }else{links_A$color[r] = paste0(paletteString, paletteTot)}}

  # to group the less "connected" links and assign them the same color (which is not the clearest one and it will be visible)
  if(count < max){num = ((paletteTot/max)*count)%/%1+1
    if(num > 2){links_A$color[r] = paste0(paletteString, num, ",z=", num)
    }else{links_A$color[r] = paste0(paletteString, "2,z=2")}
  }else{links_A$color[r] = paste0(paletteString, paletteTot, ",z=", paletteTot)}}

  # # to convert clear color coefficients to high color and viceversa (suitable only for scales from 1 to 11)
  # if(count < max){num = ((paletteTot/max)*count)%/%1+1
  #   links_A$color[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))
  # }else{num = paletteTot
  #   links_A$color[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))}}
  
  # # to regulate thickness according to level of link "connection"
  #links_A$thickness[r] = paste0("thickness=", 1+(count-1)/10, "p")}

l = dim(combinations_E)[1]
links_E = data.frame(cytoband1=rep(".",l), start1=rep(0,l), end1=rep(0,l), cytoband2=rep(".",l), start2=rep(0,l), end2=rep(0,l), color=rep(0,l), stringsAsFactors = F)
for(r in 1:l){
  gene1 = combinations_E$a[r]
  gene2 = combinations_E$b[r]
  count = combinations_E$V1[r]
  links_E$cytoband1[r] = table2[table2$symbol==gene1, "bandSTART"]
  links_E$start1[r] = table2[table2$symbol==gene1, "start_rel"]
  links_E$end1[r] = table2[table2$symbol==gene1, "end_rel"]
  links_E$cytoband2[r] = table2[table2$symbol==gene2, "bandSTART"]
  links_E$start2[r] = table2[table2$symbol==gene2, "start_rel"]
  links_E$end2[r] = table2[table2$symbol==gene2, "end_rel"]

  # # assign color according to level of link "connection"
  # if(count < max){links_E$color[r] = paste0(paletteString, ((paletteTot/max)*count)%/%1+1)
  # }else{links_E$color[r] = paste0(paletteString, paletteTot)}}

  # to group the less connected connections and assign them the same color (which is not the clearest one and it will be visible)
  if(count < max){num = ((paletteTot/max)*count)%/%1+1
    if(num > 2){links_E$color[r] = paste0(paletteString, num, ",z=", num)
    }else{links_E$color[r] = paste0(paletteString, "2,z=2")}
  }else{links_E$color[r] = paste0(paletteString, paletteTot, ",z=", paletteTot)}}

  # # to convert clear color coefficients to high color and viceversa (suitable only for scales from 1 to 11)
  # if(count < max){num = ((paletteTot/max)*count)%/%1+1
  #   links_E$color[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))
  # }else{num = paletteTot
  #   links_E$color[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))}}

  # # to regulate thickness according to level of link "connection"
  #links_E$thickness[r] = paste0("thickness=", 1+(count-1)/10, "p")}


write.table(band_circos, file=paste0(aneuplDIR, "/cytobands.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histGDcase, file=paste0(aneuplDIR, "/histGainDupAneupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histLDcase, file=paste0(aneuplDIR, "/histLossDelAneupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histGDctrl, file=paste0(aneuplDIR, "/histGainDupEupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histLDctrl, file=paste0(aneuplDIR, "/histLossDelEupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(links_A, file=paste0(aneuplDIR, "/linksAneupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(links_E, file=paste0(aneuplDIR, "/linksEupl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)




#-------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------- CHROMOTHRIPSIS ----------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# load and edit sample table
sampleID = read.delim(paste0(chromotDIR, "/arrayNames.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
selCase = sampleID[sampleID$chromothripsis==1, "Sample"] # get case samples
selCtrl = sampleID[sampleID$chromothripsis==2, "Sample"] # get ctrl samples
totCase = length(selCase) # get number of cases
totCtrl = length(selCtrl) # get number of ctrls


####


# load tables for CNA events
gainDup = read.delim(paste0(chromotDIR, "/matrixPatientGene_chromot_array_chromothripsis_gain-dup_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
lossDel = read.delim(paste0(chromotDIR, "/matrixPatientGene_chromot_array_chromothripsis_loss-del_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
# get patients list
gdsamples = gainDup$Sample
ldsamples = lossDel$Sample
# change to 1 all values higher than 1
gainDup = as.data.frame(apply(gainDup[-1], 2, function(x) sapply(x, function(y) if(y>1){y=1}else{y=y})))
lossDel = as.data.frame(apply(lossDel[-1], 2, function(x) sapply(x, function(y) if(y>1){y=1}else{y=y})))
# assign patient names to rownames
rownames(gainDup) = gdsamples
rownames(lossDel) = ldsamples

# separate samples into cases and ctrls for PATIENTxGENE matrices
gdCase = gainDup[selCase,]
gdCtrl = gainDup[selCtrl,]
ldCase = lossDel[selCase,]
ldCtrl = lossDel[selCtrl,]

# load genes list for chromothripsis
genes = read.delim(paste0(chromotDIR, "/170428_circos_genes.txt"), header=T, sep="\t", stringsAsFactors = F, check.names = F)

# annotate coordinates
genes = getCoordinatesForSymbols(genes)
# adjust end position of all genes to make them visible on plot (start + 1000000)
genes$end2 = genes$start + 500000
# edit chr names to CIRCOS format
genes$chr = paste0("hs", genes$chr)

# count events of gain+dup and loss+del for cases and ctrl and
genes$gdCase = rep(NA, dim(genes)[1])
genes$ldCase = rep(NA, dim(genes)[1])
genes$gdCtrl = rep(NA, dim(genes)[1])
genes$ldCtrl = rep(NA, dim(genes)[1])
for(r in 1:dim(genes)[1]){
  gene = genes$symbol[r]
  print(sum(ldCase[,gene]))
  if(gene %in% colnames(gdCase)){genes$gdCase[r] = sum(gdCase[,gene])/totCase*100}
  if(gene %in% colnames(ldCase)){genes$ldCase[r] = sum(ldCase[,gene])/totCase*100}
  if(gene %in% colnames(gdCtrl)){genes$gdCtrl[r] = sum(gdCtrl[,gene])/totCtrl*100}
  if(gene %in% colnames(ldCtrl)){genes$ldCtrl[r] = sum(ldCtrl[,gene])/totCtrl*100}}

# create table with gene positions for histogram track ("hist.txt")
histGDcase = genes[,c(3,4,7,8)]
histLDcase = genes[,c(3,4,7,9)]
histGDctrl = genes[,c(3,4,7,10)]
histLDctrl = genes[,c(3,4,7,11)]

# create table with gene LABELS for text track ("text.txt")
textALL = genes[,c(3,4,7,1)]


####


# load tables for segments
dgv_filtered = read.delim(paste0(chromotDIR, "/dgv_filtered.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
dgv_filtered$Sample<-unlist(lapply(strsplit(dgv_filtered$Sample, " "), "[[", 1))

# remove undesired samples
segments = dgv_filtered[which(dgv_filtered$Sample %in% sampleID$Sample),]

# remove LOH events and change other event names
segments = segments[segments$Event!="LOH",]
segments[segments$Event=="CN Loss", "Event"] = "loss"
segments[segments$Event=="CN Gain", "Event"] = "gain"
segments[segments$Event=="High Copy Gain", "Event"] = "dup"
segments[segments$Event=="Homozygous Copy Loss", "Event"] = "del"

# create event-color palette
palette = data.frame(event = c("loss", "dup",  "gain", "del"), color = c("orange", "blue", "green", "black"), stringsAsFactors = F)
palette = data.frame(event = c("loss", "dup",  "gain", "del"), color = c("lum70chr15", "lum70chr15", "lum70chr6", "lum70chr6"), stringsAsFactors = F) # for prof circos

segmentsCase = segments[segments$Sample %in% selCase$Sample,]
segmentsCtrl = segments[segments$Sample %in% selCtrl$Sample,]

# get unique samples list and calculate radius coordinates (r0 and r1) for each samples
sampleCase = unique(segmentsCase$Sample)
orderCase = data.frame(Sample=sampleCase, r0=rep(0, length(sampleCase)), r1=rep(0, length(sampleCase)), stringsAsFactors = F)
x = 0.98 # this is the outer radial limit
w = 0.0045 # this is the width of a single stripe
s = 0.004 # this is the space between the stripes
for(r in 1:length(sampleCase)){
  orderCase$r0[r] = x - w*(r-1) - s*(r-1)
  orderCase$r1[r] = x - w*(r) - s*(r-1)}

sampleCtrl = unique(segmentsCtrl$Sample)
orderCtrl = data.frame(Sample=sampleCtrl, r0=rep(0, length(sampleCtrl)), r1=rep(0, length(sampleCtrl)), stringsAsFactors = F)
x = 0.996 # this is the outer radial limit
w = 0.00165 # this is the width of a single stripe
s = 0.0003 # this is the space between the stripes
for(r in 1:length(sampleCtrl)){
  orderCtrl$r0[r] = x - w*(r-1) - s*(r-1)
  orderCtrl$r1[r] = x - w*(r) - s*(r-1)}

sampleAll = unique(segments$Sample)
orderAll = data.frame(Sample=sampleAll, r0=rep(0, length(sampleAll)), r1=rep(0, length(sampleAll)), stringsAsFactors = F)
x = 0.996 # this is the outer radial limit
w = 0.00165 # this is the width of a single stripe
s = 0.0003 # this is the space between the stripes
for(r in 1:length(sampleAll)){
  orderAll$r0[r] = x - w*(r-1) - s*(r-1)
  orderAll$r1[r] = x - w*(r) - s*(r-1)}

# extract genomic coordinate data and build datasets for HIGHLIGHTs(dataframe) and LINE PLOTs(GRanges)
case = createSegmentsForHighlight(segmentsCase, orderCase, palette)
ctrl = createSegmentsForHighlight(segmentsCtrl, orderCtrl, palette)
all = createSegmentsForHighlight(segments, orderAll, palette)

# get datasets for HIGHLIGHTs
highlightCase = case[[1]]
highlightCtrl = ctrl[[1]]
highlightAll = all[[1]]

## get datasets for LINE PLOTs
# convert GRanges to Rle (which is equivalent to calculate how many patients have a specific event in each position)
rleGainDupCase = coverage(case[[2]])
rleLossDelCase = coverage(case[[3]])
rleGainDupCtrl = coverage(ctrl[[2]])
rleLossDelCtrl = coverage(ctrl[[3]])

# subtract Rle of lossDel from Rle of gainDup and calculate percentage (on the total of patients in the group)
rleCase = ((rleGainDupCase - rleLossDelCase)/totCase)*100
rleCtrl = ((rleGainDupCtrl - rleLossDelCtrl)/totCtrl)*100

# convert back Rle to GRanges
grCase = as(rleCase, "GRanges")
grCtrl = as(rleCtrl, "GRanges")

# extract genomic coordinates and build dataset for LINE PLOT
lineCase = data.frame(chr=as.vector(seqnames(grCase)), start=start(ranges(grCase)), end=end(ranges(grCase)),  score=score(grCase))
lineCtrl = data.frame(chr=as.vector(seqnames(grCtrl)), start=start(ranges(grCtrl)), end=end(ranges(grCtrl)),  score=score(grCtrl))

write.table(histGDcase, file=paste0(chromotDIR, "/histGDCase.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histLDcase, file=paste0(chromotDIR, "/histLDCase.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histGDctrl, file=paste0(chromotDIR, "/histGDCtrl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(histLDctrl, file=paste0(chromotDIR, "/histLDCtrl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(textALL, file=paste0(chromotDIR, "/textALL.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)

write.table(highlightCase, file=paste0(chromotDIR, "/segmentsCase.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(highlightCtrl, file=paste0(chromotDIR, "/segmentsCtrl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(highlightAll, file=paste0(chromotDIR, "/segmentsAll.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)

write.table(lineCase, file=paste0(chromotDIR, "/lineCase.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)
write.table(lineCtrl, file=paste0(chromotDIR, "/lineCtrl.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)

write.table(highlightCase, file=paste0(profDIR, "/segmentsCaseProf.txt"), sep=" ", dec=".", col.names=F, row.names=F, quote=F)


#--------------------------------------------- TEST -------------------------------------------

gainDup[,"AFDN"]

gd = genes[genes$event=="gain/dup",]
ld = genes[genes$event=="loss/del",]

table(gd$symbol %in% ld$symbol)
table(ld$symbol %in% gd$symbol)


ir = IRanges(start=segmentsCase$start, end=segmentsCase$end)
resize(ranges, width=100, fix="end")
coverage(ir)

gr = GRanges(seqnames = segmentsCase$chr, ranges = ir)
coverage(gr)

gain = c(1,1,1,5,5,2,5,3,2)
loss = c(2,2,3,3,4,4)

rlGain = Rle(gain)
rlLoss = Rle(loss)

((rlGain - rlLoss)/10)*100

(11/17)*17
((11/17)*17)%/%1+1

i = 4
i+10*(1-0.2*(i-1))

summary(combinations_A$V1)
summary(combinations_E$V1)
barplot(table(combinations_A$V1))
barplot(table(combinations_E$V1))

a = c(1, 1, 1, 2, 2, 3, 4, 4)
b = c(3.5, 3.5, 2.5, 2, 2, 1, 2.2, 7)
df1 <-data.frame(a,b)
ddply(df1,.(a,b),nrow)

g = c("x","y","z","z")
df = as.data.frame(t(combn(g, 2)), stringsAsFactors=F)
colnames(df) = c("a","b")
ddply(df,.(df[,1], df[,2]),nrow)
ddply()
# grepl(band[1,2],band_circos$ID)
# 
# b = band[(band[,2]=="q24.21") & (band[,1]==12),]
# t = table2[table2$bandSTART=="12q24.21",]
# class(t)
# class(table2$start[1])
# table2$start = as.numeric(table2$start)
# class(band$V6)
# for(r in 1:dim(table2)[1]){print(class(table2$start[r]))}
# 
# band[band[,2]=="p11.23",]
# 
# # write to file
# write.table(table2, file=paste0(aneuplDIR, "/genes_data.txt"), sep="\t", dec=".", col.names=F, row.names=F, quote=F)

