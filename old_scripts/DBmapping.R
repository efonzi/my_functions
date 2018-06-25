# Copyright (C) 2015  Marco Manfrini, PhD

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#----------------------------------------------------------------------
#----------------------------------------------------------------------

# PURPOSE:  This R code is designed to annotate genes with most relevant DB (GO, Reactome, KEGG, PubMed)
#           It accepts in input a one or two columns file which contain at least gene's names. 
#           Sample column is needed to calculate the matrix of per-patients mutated pathways.
#           Columns have to be named 'SAMPLE'  and 'GENE'.

#----------------------------------------------------------------------
#----------------------------------------------------------------------

library(org.Hs.eg.db)
library(reactome.db)
library(GO.db)
#library(KEGGREST)
library(KEGG.db)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

#USAGE EXAMPLES

# columns(org.Hs.eg.db)
# 
# keytypes(org.Hs.eg.db)
# 
# mapIds(org.Hs.eg.db, keys=c("TP53"), column="PATH", keytype="SYMBOL") 
# 
# mapIds(org.Hs.eg.db, keys=c("RUNX1","TP53"), column="GO", keytype="SYMBOL", multiVals="list") 
# 
# keys(org.Hs.eg.db, keytype="SYMBOL", pattern="COX") 
# 
# select(org.Hs.eg.db, keys=c("TP53","RUNX1"), columns=c("GO","ENTREZID"), keytype="SYMBOL") 
# 
# select(GO.db, keys=c("GO:0071560"), columns=c("TERM", "ONTOLOGY", "DEFINITION"), keytype="GOID")

# keggList("pathway", "hsa") ## returns the list of human pathways

#----------------------------------------------------------------------
#----------------------------------------------------------------------

# FUNCTIONS

selAnnot <- function (GOList, ONTOLOGY) {
  sel<-which(GOList[, 'ONTOLOGY']==ONTOLOGY)
  GOList<-GOList[sel,]
  #dup<-duplicated(GOList[,1])
  #GOList<-GOList[!dup,]
  
  geneNames<-as.character(mapIds(org.Hs.eg.db, keys=GOList[,1], column="SYMBOL", keytype="ENTREZID"))
  
  #PUBMEDList<-select(org.Hs.eg.db, keys=geneNames, columns=c("PMID"), keytype="SYMBOL")
  
  GOAnno<-select(GO.db, keys=GOList[,2], column=c("TERM", "ONTOLOGY", "DEFINITION"), keytype="GOID")
  
  
  
  geneAnno<-cbind(geneNames, GOList[,1], GOList[,3], GOAnno)
  cnames<-colnames(geneAnno)
  cnames[2]<-"ENTREZID"
  cnames[3]<-"EVIDENCE"
  colnames(geneAnno)<-cnames
  return(geneAnno)
}

selAnnotKEGG <- function (KEGGList) {
  sel<-which(!is.na(KEGGList$PATH))
  KEGGList<-KEGGList[sel,]
  dup<-duplicated(KEGGList[,1])
  KEGGList<-KEGGList[!dup,]
  
  KEGGpaths<-unlist(as.list(KEGGPATHID2NAME))
  KEGGpathCol<-rep("", dim(KEGGList)[1])
  
  for(i in 1:dim(KEGGList)[1]){
    sel<-which(names(KEGGpaths)==KEGGList[i, 'PATH'])
    KEGGpathCol[i]<-KEGGpaths[sel]
  }
  
  KEGGList<-cbind(KEGGList, KEGGpathCol)
  
  colnames(KEGGList)<-c("GENE", "PAHID", "PATHNAME")
  
  return(KEGGList)
}


makePubMedTable<-function(PUBMEDList){
  
  items<-levels(as.factor(PUBMEDList$SYMBOL))
  
  PUBMEDListTab=NULL
    
  for(i in 1:length(items)){
    
    PUBMEDIds<-paste(PUBMEDList[PUBMEDList$SYMBOL==items[i], 'PMID'], collapse=", ")
    geneSymbol<-items[i]
    newRow<-cbind(geneSymbol, PUBMEDIds)
    PUBMEDListTab<-rbind(PUBMEDListTab, newRow)
    
  }
  
  return(PUBMEDListTab)
  
}

makePatKEGGPathsMatrix<-function(perSampleMutFileName){
    
  sel<-which(!is.na(KEGGList$PATH))
  KEGGList<-KEGGList[sel,]
  dup<-duplicated(KEGGList[,1])
  KEGGList<-KEGGList[!dup,]
  dup<-duplicated(KEGGList[,2])
  KEGGList<-KEGGList[!dup,]
  
  head(perSampleMutFileName)
  perSampleMut<-read.delim(paste(getwd(), perSampleMutFileName, sep="/"), header=T, stringsAsFactors=F)
  nsamples<-length(levels(as.factor(perSampleMut$SAMPLE)))
  IDsamples<-levels(as.factor(perSampleMut$SAMPLE))
  
  ncol<-dim(KEGGList)[1]
  nrow<-nsamples
  
  m<-matrix(rep(0, ncol*nrow), ncol=ncol, byrow=T)
  
  tmp<-data.frame()
  
  for(i in 1:nrow){
    
    patsID<-IDsamples[i]
    tmp<-perSampleMut[perSampleMut$SAMPLE==patsID,]
    for(j in 1:dim(tmp)[1]){
      
      gene<-tmp[j, 'GENE']
      sel<-which(KEGGList$SYMBOL==gene)
      
      if(length(sel)>0){
      
        m[i, sel]<-1
      }
      
    }
      
  }
  
  colnames(m)<-KEGGList$PATH
  rownames(m)<-IDsamples
  
  return(m)
  
}


#----------------------------------------------------------------------
#----------------------------------------------------------------------

#setwd("/home/marco/Dati/Bologna/Anna/dbmapping")
#setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/1611_tp53")

#samplesFileName="genes_list.tsv"   # THIS INPUT FILES SHOULD BE A 2 COLUMN TABLE WITH 'SAMPLE' AND 'GENE' INFO
                              # FOR ANNOTATION ONLY IT IS ENOUGH A 1 COLUMN FILE WITH 'GENE' NAMES
#dset<-read.delim(paste(getwd(), samplesFileName, sep="/"), header=T, stringsAsFactors=F)

setwd(dirname(commandArgs()[6]))
samplesFileName=basename(commandArgs()[6])
dset<-read.delim(samplesFileName, header=T, stringsAsFactors=F)

#patsName<-dset[, 'SAMPLE']

geneNames<-dset[, 'GENE']

ENTREZList<-as.character(mapIds(org.Hs.eg.db, keys=geneNames, column="ENTREZID", keytype="SYMBOL"))

ENZYMEList<-as.character(mapIds(org.Hs.eg.db, keys=geneNames, column="ENZYME", keytype="SYMBOL"))

PATHList<-as.character(mapIds(org.Hs.eg.db, keys=geneNames, column="PATH", keytype="SYMBOL"))

PUBMEDList<-select(org.Hs.eg.db, keys=geneNames, column="PMID", keytype="SYMBOL")

REACTOMEList<-select(reactome.db, keys=ENTREZList, columns=c("PATHNAME", "PATHID"), keytype=c("ENTREZID"))

GOList<-select(org.Hs.eg.db, keys=ENTREZList, columns=c("GO"), keytype="ENTREZID")

KEGGList<-select(org.Hs.eg.db, keys=geneNames, columns=c("PATH"), keytype="SYMBOL")

geneAnnoBP<-selAnnot(GOList, "BP")
geneAnnoMF<-selAnnot(GOList, "MF")
geneAnnoCC<-selAnnot(GOList, "CC")

geneAnnoPubMed<-makePubMedTable(PUBMEDList)

geneAnnoKEGG<-selAnnotKEGG(KEGGList)

##########################################################
#pats2KEGGMatrix<-makePatKEGGPathsMatrix(samplesFileName)
##########################################################

#write.table(ENTREZList, file="ENTREZList.tsv", sep="\t", row.names=F, col.names=T, quote=F)
write.table(ENTREZList, file=paste("ENTREZList", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)

write.table(geneAnnoBP, file=paste("geneAnnoBP", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(geneAnnoMF, file=paste("geneAnnoMF", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(geneAnnoCC, file=paste("geneAnnoCC", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(geneAnnoPubMed, file=paste("geneAnnoPubMed", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(geneAnnoKEGG, file=paste("geneAnnoKEGG", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
#write.table(pats2KEGGMatrix, file=paste("pats2KEGGMatrix", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=T, col.names=NA, quote=F)

write.table(ENZYMEList, file=paste("ENZYMEList", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(PATHList, file=paste("PATHList", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(PUBMEDList, file=paste("PUBMEDList", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
write.table(REACTOMEList, file=paste("REACTOMEList", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)

#write.table(ENZYMEList, file="ENZYMEList.tsv", sep="\t", row.names=F, col.names=T, quote=F)
#write.table(PATHList, file="PATHList.tsv", sep="\t", row.names=F, col.names=T, quote=F)
#write.table(PUBMEDList, file="PUBMEDList.tsv", sep="\t", row.names=F, col.names=T, quote=F)
#write.table(REACTOMEList, file="REACTOMEList.tsv", sep="\t", row.names=F, col.names=T, quote=F)

