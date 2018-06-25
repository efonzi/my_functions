#-------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------- Create pathID-pathName-entrezID-symbol table for REACTOME and KEGG ---------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# launch as R script
# first argument must be the working directory
# second argument must be the name of the reactome hierarchy file to be used

library(org.Hs.eg.db)
library(reactome.db)
library(KEGGREST)

pathDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/pathway_databases"

cat("retrieving all pathway and gene info for REACTOME\n")
# load and edit REACTOME hierarchy file
reac_hierarchy = read.delim(paste0(pathDIR, "/reactome_hierarchy_human.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
reac_hierarchy = as.data.frame(sapply(reac_hierarchy, as.character), stringsAsFactors=F)
# get all REACTOME pathIDs
reactome_path = keys(reactome.db, keytype="PATHID")
# create pathID-to-entrezID table
reactome_path = select(reactome.db, keys=reactome_path, column="ENTREZID", keytype = "PATHID")
# get all pathNames from pathIDs
pathNames = select(reactome.db, keys=reactome_path$PATHID, column="PATHNAME", keytype = "PATHID")[,2]
# get all symbols from entrezIDs
symbols = select(org.Hs.eg.db, keys=reactome_path$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
# make pathID-pathName-entrezID-symbol table
reactome_path = data.frame(pathID=reactome_path$PATHID, pathName=pathNames, entrezID=reactome_path$ENTREZID, symbol=symbols, stringsAsFactors=F)
# select only human pathways
reactome_path = reactome_path[which(grepl("Homo sapiens", reactome_path$pathName)),]
# remove "Homo sapiens: " prefix from pathNames
reactome_path$pathName = substr(reactome_path$pathName, 15, nchar(reactome_path$pathName))

write.table(reactome_path, file=paste0(pathDIR, "/pathID_pathName_entrezID_symbol_REACTOME.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)

cat("retrieving all pathway and gene info for KEGG\n")
# load KEGG hierarchy file
kegg_hierarchy = read.delim(paste0(pathDIR, "/kegg_pathway_hierarchy_edited.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
# get all human ("hsa") KEGG pathIDs (removing "path" prefix)
kegg_path = sapply(names(keggList("pathway", "hsa")), function(x) strsplit(x, ":")[[1]][2])
# get entrezIDs and symbols for each pathway
genes = sapply(kegg_path, function(x) keggGet(x)[[1]]$GENE) # it's a list of lists (each has, for a pathway, entrezids at odd indices - symbols at even indices)
# extract entrezIDs (obtains another list of lists)
entrezID = sapply(genes, function(x) if(!is.null(x)){x[seq(1,length(x),2)]}) # for each list, select only odd indices (entrez ids)
names(entrezID) = sapply(names(entrezID), function(x) strsplit(x, ":")[[1]][2]) # remove "path:" prefix from pathIDs
# extract symbols (obtains another list of lists)
symbols = sapply(genes, function(x) if(!is.null(x)){x[seq(2,length(x),2)]}) # for each list, select only even indices (symbols)
symbols = sapply(symbols, function(y) as.vector(sapply(y, function(x) strsplit(x, ";")[[1]][1]))) # get symbols (prefix of a longer string)
names(symbols) = sapply(names(symbols), function(x) strsplit(x, ":")[[1]][2]) # remove "path:" prefix from pathIDs
# make pathID-pathName-entrezID-symbol table
kegg_path=NULL
for(l in 1:length(symbols)){
  pathID = names(symbols)[l]
  pathName = kegg_hierarchy[kegg_hierarchy$pathID==names(symbols)[l], "pathName"]
  gene_number = length(symbols[[l]])
  entr = entrezID[[l]]
  symb = symbols[[l]]         
  kegg_path = rbind(kegg_path, data.frame(pathID=rep(pathID, gene_number), pathName=rep(pathName, gene_number), entrezID=entr, symbol=symb, stringsAsFactors = F))
}

write.table(kegg_path, file=paste0(pathDIR, "/pathID_pathName_entrezID_symbol_KEGG.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------- 170421 - KEGG hierarchy file editing ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

pathDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/pathway_databases"

# load hierarchy file
hierarchy = read.delim(paste0(pathDIR, "/kegg_pathway_hierarchy.txt"), header=F, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
# remove unwanted lines and columns
hierarchy = hierarchy[9:590,]
# remove lines with "#"
hierarchy = hierarchy[grep("#", hierarchy$V1, invert=T),]

# get number of rows and add columns
l = dim(hierarchy)[1]
hierarchy = data.frame(hierarchy[,1], rep(NA,l), rep(NA, l), rep(NA, l), rep(NA, l), stringsAsFactors=F)
colnames(hierarchy) = c("level", "pathID", "pathName", "level_B", "level_A")

# remove "<b>" and "</b>" from A level names
sel = grep("A<", hierarchy$level)
hierarchy[sel, "level"] = sapply(hierarchy[sel, "level"], function(x) strsplit(x, "</b>")[[1]][1])
hierarchy[sel, "level"] = sapply(hierarchy[sel, "level"], function(x) paste(strsplit(x, "<b>")[[1]][1], strsplit(x, "<b>")[[1]][2], collapse = " "))

# select A levels and write "A" in "level" columns and name in "pathName"
selA = which(sapply(hierarchy$level, function(x) strsplit(x, " ")[[1]][1] == "A"))
hierarchy[selA, c("level", "pathName")] = c(sapply(hierarchy[selA, "level"], function(x) strsplit(x, " ")[[1]][1]), sapply(hierarchy[selA, "level"], function(x) strsplit(x, "A ")[[1]][2]))

# select B levels and write "B" in "level" columns and name in "pathName"
selB = which(sapply(hierarchy$level, function(x) strsplit(x, " ")[[1]][1] == "B"))
hierarchy[selB, c("level", "pathName")] = c(sapply(hierarchy[selB, "level"], function(x) strsplit(x, " ")[[1]][1]), sapply(hierarchy[selB, "level"], function(x) strsplit(x, "B ")[[1]][2]))

# select C levels and write "C" in "level" columns, ID in "pathID" and name in "pathName"
# (for 2nd column we need to select the 5th element of the split string [[1]][5] because there are 4 spaces between "C" and pathID)
selC = which(sapply(hierarchy$level, function(x) strsplit(x, " ")[[1]][1] == "C"))
hierarchy[selC, c("level", "pathID", "pathName")] = c(sapply(hierarchy[selC, "level"], function(x) strsplit(x, " ")[[1]][1]), sapply(hierarchy[selC, "level"], function(x) strsplit(x, " ")[[1]][5]), hierarchy[selC, "level"])
# subset pathNames (remove level and pathID info)
hierarchy[selC, "pathName"] = substr(hierarchy[selC, "pathName"], 13, nchar(hierarchy[selC, "pathName"]))

# assign level A and B names to each pathway
for(r in 1:l){
  if(hierarchy[r,"level"] == "A"){hierarchy$level_A[r:l] = hierarchy[r,"pathName"]}
  if(hierarchy[r,"level"] == "B"){hierarchy$level_B[r:l] = hierarchy[r,"pathName"]}
}

# subset for C-level pathways and remove "level" column
hierarchy = hierarchy[hierarchy$level == "C", -1]

# add "hsa" to pathIDs
hierarchy$pathID = paste0("hsa", hierarchy$pathID)

write.table(hierarchy, file=paste0(pathDIR, "/kegg_pathway_hierarchy_edited.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)


#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------- 170419 - REACTOME hierarchy file editing ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
# load hierarchy file
hierarchy = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/pathway_databases/reactome_hierarchy.tsv", header=F, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "na")
# get only human pathways
hierarchy = hierarchy[grep("HSA", hierarchy[,1]),]
# extract pathID (remove "R-HSA-")
hierarchy = as.data.frame(apply(hierarchy, 2, function(x) sapply(x, function(x) strsplit(x, "-")[[1]][3])), stringsAsFactors=F)
# change colnames
colnames(hierarchy) = c("parent", "child")

write.table(hierarchy, file="/Users/Eugenio/Desktop/Lavoro/Seragnoli/pathway_databases/reactome_hierarchy_human.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)


#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------ 170411 - Brewer color palette test for Aneuploidy CIRCOS legend ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)

pal = brewer.pal(9, "Blues")

### to use Color Brewer as discrete colors
# create dummy dataset and save color numbers as character
# as we are using colors as discrete entities, we must calculate the exact number for each color in the palette
data = data.frame(one=1:17, two=1:17)
data$col = as.character(((9/dim(data)[1])*data[,1])%/%1+1)
data$col[dim(data)[1]] = as.character(9)
ggplot(data, aes(x=one, y=two, color=col)) + geom_point() + scale_color_brewer(palette="Blues")

### to use Color Brewer as continuous colors
# create dummy dataset and save color numbers as integer/numeric
# as we are using colors as continuous entities, we can give a level for each data point OR
data = data.frame(one=1:17, two=1:17, col=1:17)
ggplot(data, aes(x=one, y=two, color=col)) + geom_point() + scale_color_gradientn(colours=pal)
# we can calculate the exact number for each color in the palette
data$col = ((9/dim(data)[1])*data[,1])%/%1+1
data$col[dim(data)[1]] = 9
ggplot(data, aes(x=one, y=two, color=col)) + geom_point() + scale_color_gradientn(colours=pal)



#-------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------- 1703 - mutations and CNA on MCL1 ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
snpIndel = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis/wesTable_SNVindel.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "na")
wesID_snpID = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/wesID_snpID.txt", header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "na")
sampleID = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis/arrayNames.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "na")
sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
wesID_snpID$snpID<-unlist(lapply(strsplit(wesID_snpID$snpID, " "), "[[", 1))

g = c("MYC", "BTRC", "MYCBP2", "PLK1", "PTTG1", "BUB1", "MAD2L1", "BUB1B", "HUWE1", "TRIM17", "USP9X", "MCL1", "SKP1", "CUL1", "FBXW7", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC16", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "CDC16", "CDC20B", "CDC23", "CDC27")

match=NULL
for(r in 1:dim(snpIndel)[1]){
  print(r)
  if(snpIndel[r,"GENE"] %in% g){
    wesID = snpIndel[r,"Sample"]
    snpID = wesID_snpID[wesID_snpID$wesID == wesID, "snpID"]
    if(!is.na(snpID)){
      chromotripsis = sampleID[sampleID$Sample == snpID, "chromothripsis"]
      if(chromotripsis == 1){chromotripsis = "case"}
      if(chromotripsis == 2){chromotripsis = "ctrl"}
    }else{chromotripsis = NA}
    df = data.frame(snpIndel[r,c(1,6)], snpID=snpID, chromotripsis=chromotripsis)
    match = rbind(match, df)}}

write.table(match, file="/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/mcl1&co.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
