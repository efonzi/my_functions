
library(org.Hs.eg.db)
library(reactome.db)
library(clusterProfiler)
library(ReactomePA)
library(GO.db)
library(KEGGREST)




# loop through each patient and perform hypergeometric test
for(r in 1:dim(geneTable)[1]){
  
  genes = colnames(geneTable)[which(geneTable[1,] == 1)]
  # genes = genes[complete.cases(genes)] <---- I shouldn't have NA's in the genes list, so I don't need this step
  ids = as.vector(unlist(mget(genes,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
  ids = ids[complete.cases(ids)]

  ## reactome
  reactome = enrichPathway(ids, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1)
  reactome = as.data.frame(reactome)
  
  kegg = enrichKEGG(ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  
  go = enrichGO(gene=ids, OrgDb=org.Hs.eg.db, keytype = "ENTREZID", ont = "BP", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1)
  
}

as.data.frame(kegg)
as.data.frame(go)




r = as.vector(unlist(lapply(reactome@geneSets, length)))
k = as.vector(unlist(lapply(kegg@geneSets, length)))
g = as.vector(unlist(lapply(go@geneSets, length)))
hist(log10(r), breaks=length(r))
hist(log10(k), breaks=length(k))
hist(log10(g), breaks=length(g))


#### REACTOME

# GET ALL REACTOME PATHIDs



# universe reactome
univ=NULL
for(i in 1:length(reactome@geneSets)){univ = unique(append(univ, reactome@geneSets[[i]]))}
reactome@universe[which(!reactome@universe %in% univ)]

# pathID reactome from reactome PA
pathID=names(reactome@geneSets)

reac = keys(reactome.db, keytype = "PATHNAME")
reac[grep("Homo", reac)]
reacID = keys(reactome.db, keytype = "PATHID")

table(pathID %in% reacID)
table(reacID %in% pathID)
reacIDno = reacID[which(!reacID %in% pathID)]
reacNameNo = mapIds(reactome.db, keys=reacIDno, column="PATHNAME", keytype = "PATHID")
reacNameNo = reacNameNo[grep("Homo", reacNameNo)]

pathName = mapIds(reactome.db, keys=pathID, column="PATHNAME", keytype = "PATHID")
pathName[which(duplicated(pathName))]
pathName[grep("RNA Pol II", pathName)]

#### KEGG
pathID=names(kegg@geneSets)

kegg = keys(org.Hs.eg.db, keytype = "PATH")

#reacID = keys(reactome.db, keytype = "PATHID")
columns(org.Hs.eg.db)
############
KEGGPATHID2NAME
KEGGpaths<-unlist(as.list(KEGGPATHID2NAME))
KEGGpathCol<-rep("", dim(KEGGList)[1])
k = as.list(KEGGPATHID2NAME)

###############

#### GOBP

#----------------------------------------------------------------------------------------------------------------------------------------------

# genes = colnames(gain)[which(gain[3,] == 1)]
# ids = as.vector(unlist(mget(genes,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
# ids = ids[complete.cases(ids)]
# reactome = enrichPathway(ids, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1)
# kegg = enrichKEGG(ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
# go = enrichGO(gene=ids, OrgDb=org.Hs.eg.db, keytype = "ENTREZID", ont = "BP", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1)
# rSets = names(reactome@geneSets)
# 
# 
# genes = colnames(gain)[which(gain[4,] == 1)]
# ids = as.vector(unlist(mget(genes,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
# ids = ids[complete.cases(ids)]
# reactome = enrichPathway(ids, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1)
# kegg = enrichKEGG(ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
# go = enrichGO(gene=ids, OrgDb=org.Hs.eg.db, keytype = "ENTREZID", ont = "BP", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1)


setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/Dropbox/RNAseq_AML_ALL/ALL_SNP/170225_arrayTables")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy")
system(paste0("ls ",getwd()))
del = read.delim("arrayTable_del.tsv", stringsAsFactors = F, sep = "\t")
dup = read.delim("arrayTable_dup.tsv", stringsAsFactors = F, sep = "\t")
loss = read.delim("arrayTable_loss.tsv", stringsAsFactors = F, sep = "\t")
gain = read.delim("arrayTable_gain.tsv", stringsAsFactors = F, sep = "\t")

samples = read.delim("arrayNames.tsv", stringsAsFactors = F, sep = "\t")
