# Copyright (C) 2015-2016  Marco Manfrini, PhD

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

#TODOLIST:  CREARE LA FUNCTION PER L'ANNOTAZIONE. AL MOMENTO E' NECESSARIO ESEGUIRE IL CODICE UNA VOLTA PER OGNI GRUPPO DI GENI CAMBIANDO 
#           MANUALMENTE IL NOME DEI FILES DI OUTPUT.

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

makePatKEGGPathsMatrix<-function(parSampleMutFileName){
    
  sel<-which(!is.na(KEGGList$PATH))
  KEGGList<-KEGGList[sel,]
  dup<-duplicated(KEGGList[,1])
  KEGGList<-KEGGList[!dup,]
  dup<-duplicated(KEGGList[,2])
  KEGGList<-KEGGList[!dup,]
  
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

findsig<-function(set, alpha, outfile, nchr, trim){
  modup<-glm(Freq~Var1, data=set, family=poisson(link = "log"))
  pvals<-coef(summary(modup))[,4]
  pvals<-pvals[2:length(pvals)]
  qvals<-p.adjust(pvals, method="BH", n=length(pvals))
  p_names=substring(names(pvals), nchr, nchar(names(pvals)))
  sigs_paths<-data.frame(p_names, pvals, qvals)
  sigs_paths<-sigs_paths[order(qvals),]
  
  sel <- which(sigs_paths$qvals<alpha)
  p_names <- sigs_paths[sel, 'p_names']
  if(trim>0){
    set$Var1<-substring(set$Var1, trim, nchar(as.character(set$Var1)))
  }
  write.table(sigs_paths[sel,], file=paste0(outfile, "_0.05.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(set, file=paste0(outfile, "_freq.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  
}

toreport<-function(index){
  #write.table(ENTREZList, file=paste0("ENTREZList", "_", index, ".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(geneAnnoBP, file=paste0("geneAnnoBP", "_", index, ".tsv"), sep="\t", row.names=F, col.names=F, quote=F)
  #write.table(geneAnnoMF, file=paste0("geneAnnoMF", "_", index, ".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  #write.table(geneAnnoCC, file=paste0("geneAnnoCC", "_", index, ".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  #write.table(geneAnnoPubMed, file=paste("geneAnnoPubMed", index, unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(geneAnnoKEGG, file=paste0("geneAnnoKEGG", "_", index, ".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(REACTOMEList, file=paste0("geneAnnoReactome", "_", index, ".tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  #write.table(pats2KEGGMatrix, file=paste("pats2KEGGMatrix_A", unlist(strsplit(samplesFileName, split="[.]"))[1], "tsv", sep="."), sep="\t", row.names=T, col.names=NA, quote=F)
}

# Eugenio functions
fillPathTable <- function(allPathTable, alteredPathTable, patient){
  
  for(i in 1:length(alteredPathTable)){
    allPathTable[patient, alteredPathTable[i]] = 1
  }
  return(allPathTable)
  
}

pathwayCaseCtrlComparison <- function(pathwayTable, selCase, selCtrl){
  
  pathwayCase = pathwayTable[which(pathwayTable$Samples %in% selCase$CHIPID),]
  pathwayCtrl = pathwayTable[which(pathwayTable$Samples %in% selCtrl$CHIPID),]
  
  caseCounts<-apply(pathwayCase[,2:dim(pathwayCase)[2]], 2, sum)
  ctrlCounts<-apply(pathwayCtrl[,2:dim(pathwayCtrl)[2]], 2, sum)
  
  Fmatrix<-data.frame(caseCounts, ctrlCounts)
  sumCase<-dim(pathwayCase)[1]
  sumCtrl<-dim(pathwayCtrl)[1]
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATHWAY", "%Case", "%Ctrl", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    cat(i, "\n")
    posCase<-Fmatrix[i, 'caseCounts']
    posCtrl<-Fmatrix[i, 'ctrlCounts']
    negCase<-sumCase-posCase
    negCtrl<-sumCtrl-posCtrl
    
    x<-matrix(c(posCase, posCtrl, negCase, negCtrl), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 2]<-posCase/sumCase*100
      Rmatrix[i, 3]<-posCtrl/sumCtrl*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$PATHWAY<-rownames(Fmatrix)
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=5)
  return(Rmatrix)
  
}


#--------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------ PATHWAY ANNOTATION AND ENRICHMENT - FOR CNV/WES -------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

#setwd("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/tables_revised")
#setwd("/home/mmanfrini/data1/data/Bologna/Chromotripsy/last")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy") # FOR ARRAY FOLDER

geneCounts = list.files(getwd(), "Genes")

for(i in 1:length(geneCounts)){
  print(i) ##################################
  samplesFileName = geneCounts[i]   # THIS INPUT FILES SHOULD BE A 2 COLUMN TABLE WITH 'GENE' AND 'COUNTS' INFO
  
  index=paste0(strsplit(samplesFileName, "G")[[1]][1], "_", substr(samplesFileName, nchar(samplesFileName)-4, nchar(samplesFileName)-4))
  
  dset<-read.delim(paste(getwd(), samplesFileName, sep="/"), header=T, stringsAsFactors=F)
  colnames(dset) = c("GENE", "COUNTS")
  dset<-dset[which(dset$COUNTS!=0),]
  
  sel<-grep("LOC", dset$GENE)
  if(length(sel)>0) dset<-dset[-sel,]
  sel<-grep("LINC", dset$GENE)
  if(length(sel)>0) dset<-dset[-sel,]
  sel<-grep("MIR", dset$GENE)
  if(length(sel)>0) dset<-dset[-sel,]
  sel<-grep("SNORD", dset$GENE)
  if(length(sel)>0) dset<-dset[-sel,]
  
  geneNames<-dset[, 'GENE']
  geneNames<-gsub("[.]", "-", geneNames)
  
  # INFOS RETRIVIAL
  
  ENTREZList<-as.character(mapIds(org.Hs.eg.db, keys=geneNames, column="ENTREZID", keytype="SYMBOL"))
  REACTOMEList<-select(reactome.db, keys=ENTREZList, columns=c("PATHNAME", "PATHID"), keytype=c("ENTREZID"))
  GOList<-select(org.Hs.eg.db, keys=ENTREZList, columns=c("GO"), keytype="ENTREZID")
  KEGGList<-select(org.Hs.eg.db, keys=geneNames, columns=c("PATH"), keytype="SYMBOL")
  # ENZYMEList<-as.character(mapIds(org.Hs.eg.db, keys=geneNames, column="ENZYME", keytype="SYMBOL"))
  # PATHList<-as.character(mapIds(org.Hs.eg.db, keys=geneNames, column="PATH", keytype="SYMBOL"))
  # PUBMEDList<-select(org.Hs.eg.db, keys=geneNames, column="PMID", keytype="SYMBOL")
  
  geneAnnoBP<-selAnnot(GOList, "BP")
  geneAnnoKEGG<-selAnnotKEGG(KEGGList)
  sel<-grep("Bos taurus", REACTOMEList[,2])
  if(length(sel>0)) { REACTOMEList<-REACTOMEList[-sel,] }
  sel<-grep("Mus musculus", REACTOMEList[,2])
  if(length(sel>0)) { REACTOMEList<-REACTOMEList[-sel,] }
  # geneAnnoPubMed<-makePubMedTable(PUBMEDList)
  # geneAnnoMF<-selAnnot(GOList, "MF")
  # geneAnnoCC<-selAnnot(GOList, "CC")
  # pats2KEGGMatrix<-makePatKEGGPathsMatrix(samplesFileName)
  
  ## ENRICHMENT
  
  pset<-as.data.frame(table(REACTOMEList$PATHNAME), stringsAsFactors=F)
  gset<-as.data.frame(table(geneAnnoBP$TERM))
  kset<-as.data.frame(table(geneAnnoKEGG$PATHNAME))
  # gmfset<-as.data.frame(table(geneAnnoMF$TERM))
  # gccset<-as.data.frame(table(geneAnnoCC$TERM))
    
  findsig(set=pset, alpha=1, outfile=paste0("Reactome_", index), nchr=19, trim=15)
  findsig(set=gset, alpha=1, outfile=paste0("GOBP_", index), nchr=5, trim=0)
  findsig(set=kset, alpha=1, outfile=paste0("KEGG_", index),  nchr=5, trim=0)
  # findsig(set=gmfset, alpha=1, outfile=paste0("GOMF_", index),  nchr=5, trim=0)
  # findsig(set=gccset, alpha=1, outfile=paste0("GOCC_", index),  nchr=5, trim=0)
  
  
  ## REPORTING
  
  toreport(index)
  
}


#--------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------- PATIENT-PATHWAY ENRICHMENT ----------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy")

tablesFiles = list.files(getwd(), "AllCountsPerPatientTable_filtered")

for(i in 1:length(tablesFiles)){
  
  tableFileName = tablesFiles[i]
  index = strsplit(tableFileName, "CountsPer")[[1]][1]
  
  countTable<-read.delim(tableFileName, header=T, sep="\t", dec=".", stringsAsFactors = F)
  countTable$Sample<-unlist(lapply(strsplit(countTable$Sample, " "), "[[", 1))
  rownames(countTable) = countTable[,1]
  colnames(countTable)<-gsub("[.]", "-", colnames(countTable))
  countTable = countTable[,-1]
  
  countTableGenes = colnames(countTable)
  countTablePats = rownames(countTable)
  
  ENTREZList<-as.character(mapIds(org.Hs.eg.db, keys=countTableGenes, column="ENTREZID", keytype="SYMBOL"))
  GOList<-select(org.Hs.eg.db, keys=ENTREZList, columns=c("GO"), keytype="ENTREZID")
  REACTOMEList<-select(reactome.db, keys=ENTREZList, columns=c("PATHNAME", "PATHID"), keytype=c("ENTREZID"))
  KEGGList<-select(org.Hs.eg.db, keys=countTableGenes, columns=c("PATH"), keytype="SYMBOL")
   
  geneAnnoBP<-selAnnot(GOList, "BP")
  geneAnnoKEGG<-selAnnotKEGG(KEGGList)
  sel<-grep("Bos taurus", REACTOMEList[,2])
  if(length(sel>0)) { REACTOMEList<-REACTOMEList[-sel,] }
  sel<-grep("Mus musculus", REACTOMEList[,2])
  if(length(sel>0)) { REACTOMEList<-REACTOMEList[-sel,] }
  
  pathwayREACTOME <-sort(unique(REACTOMEList$PATHNAME))
  pathwayREACTOME = substr(pathwayREACTOME, 15, nchar(pathwayREACTOME))
  pathwayTableREACTOME <-data.frame(matrix(rep(0, length(countTablePats)*length(pathwayREACTOME)), ncol=length(pathwayREACTOME)))
  colnames(pathwayTableREACTOME) = pathwayREACTOME
  rownames(pathwayTableREACTOME) = countTablePats
  
  pathwayGOBP<-sort(unique(geneAnnoBP$TERM))
  pathwayTableGOBP <-data.frame(matrix(rep(0, length(countTablePats)*length(pathwayGOBP)), ncol=length(pathwayGOBP)))
  colnames(pathwayTableGOBP) = pathwayGOBP
  rownames(pathwayTableGOBP) = countTablePats

  pathwayKEGG<-sort(unique(as.character(geneAnnoKEGG$PATHNAME)))
  pathwayTableKEGG <-data.frame(matrix(rep(0, length(countTablePats)*length(pathwayKEGG)), ncol=length(pathwayKEGG)))
  colnames(pathwayTableKEGG) = pathwayKEGG
  rownames(pathwayTableKEGG) = countTablePats
  
  for(j in 1:length(countTablePats)){
    
    patient = countTablePats[j]
    cat(j, index, "- annotating pathways for patient", patient, "\n")
    
    sel = which(countTable[patient,] > 0)
    altered_genes = unique(colnames(countTable[sel]))
    
    ENTREZList<-as.character(mapIds(org.Hs.eg.db, keys=altered_genes, column="ENTREZID", keytype="SYMBOL"))
    GOList<-select(org.Hs.eg.db, keys=ENTREZList, columns=c("GO"), keytype="ENTREZID")
    altered_pathwayGOBP<-sort(unique(selAnnot(GOList, "BP")[,5]))
    pathwayTableGOBP = fillPathTable(pathwayTableGOBP, altered_pathwayGOBP, patient)

    KEGGList<-select(org.Hs.eg.db, keys=altered_genes, columns=c("PATH"), keytype="SYMBOL")
    altered_pathwayKEGG<-sort(unique(as.character(selAnnotKEGG(KEGGList)[,3])))
    if(length(altered_pathwayKEGG) > 0){
      pathwayTableKEGG = fillPathTable(pathwayTableKEGG, altered_pathwayKEGG, patient)}

    ENTREZList<-as.character(mapIds(org.Hs.eg.db, keys=altered_genes, column="ENTREZID", keytype="SYMBOL"))
    if("TRUE" %in% (ENTREZList %in% keys(reactome.db))){
      REACTOMEList<-select(reactome.db, keys=ENTREZList, columns=c("PATHNAME", "PATHID"), keytype=c("ENTREZID"))
      sel<-grep("Bos taurus", REACTOMEList[,2])
      if(length(sel>0)) { REACTOMEList<-REACTOMEList[-sel,] }
      sel<-grep("Mus musculus", REACTOMEList[,2])
      if(length(sel>0)) { REACTOMEList<-REACTOMEList[-sel,] }
      altered_pathwayREACTOME<-sort(unique(substr(REACTOMEList[,2], 15, nchar(REACTOMEList[,2]))))
      if(length(altered_pathwayREACTOME) > 0){
        pathwayTableREACTOME = fillPathTable(pathwayTableREACTOME, altered_pathwayREACTOME, patient)}
    }  
  }
  
  pathwayTableGOBP = cbind(Samples=countTablePats, pathwayTableGOBP)
  pathwayTableKEGG = cbind(Samples=countTablePats, pathwayTableKEGG)
  pathwayTableREACTOME = cbind(Samples=countTablePats, pathwayTableREACTOME)
  
  write.table(pathwayTableGOBP, file=paste0(index, "PathwayCountTableGOBP.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(pathwayTableKEGG, file=paste0(index, "PathwayCountTableKEGG.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(pathwayTableREACTOME, file=paste0(index, "PathwayCountTableREACTOME.tsv"), sep="\t", row.names=F, col.names=T, quote=F)
  
  chipid <- read.delim("chipnames", stringsAsFactors=FALSE)
  selCase<-chipid[which(chipid$AvsE=="1" & chipid$CHIPID!=""),]
  selCtrl<-chipid[which(chipid$AvsE=="2" & chipid$CHIPID!=""),]
  
  #execute function for fisher's test on differential pathway enrichment between case and ctrl
  GOBPcomparison = pathwayCaseCtrlComparison(pathwayTableGOBP, selCase, selCtrl)
  KEGGcomparison = pathwayCaseCtrlComparison(pathwayTableKEGG, selCase, selCtrl)
  REACTOMEcomparison = pathwayCaseCtrlComparison(pathwayTableREACTOME, selCase, selCtrl)
  
  write.table(GOBPcomparison, file=paste0(index, "_GOBP_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  write.table(KEGGcomparison, file=paste0(index, "_KEGG_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  write.table(REACTOMEcomparison, file=paste0(index, "_REACTOME_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
}


########################################################


#------------------------------------------------------------------------------------------------------------------
#-------------------------------------------- FIND EXCLUSIVE PATHWAYS BETWEEN GROUPS ------------------------------
#------------------------------------------------------------------------------------------------------------------

label = c("gainAll", "lossAll")

for(i in 1:length(label)){
  
  print(i) #################################
  index = label[i]
  
  # REACTOME
  
  rA<-read.delim(paste0(getwd(), "/Reactome_", index, "_A_0.05.tsv"), header=T, sep="\t", dec=".", stringsAsFactors=F) # AGGIUNGI LOOP X FARE TUTTI I FILE
  rfA<-read.delim(paste0(getwd(), "/Reactome_", index, "_A_freq.tsv"), header=T, sep="\t", dec=".", stringsAsFactors=F)
  
  rE<-read.delim(paste0(getwd(), "/Reactome_", index, "_E_0.05.tsv"), header=T, sep="\t", dec=".", stringsAsFactors=F)
  rfE<-read.delim(paste0(getwd(), "/Reactome_", index, "_E_freq.tsv"), header=T, sep="\t", dec=".", stringsAsFactors=F)
  
  sel<-which(rA$p_names %in% rE$p_names)
  rAUniqueList<-rA[-sel,]
  write.table(rAUniqueList, file=paste0("Reactome_unique_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  sel<-which(rE$p_names %in% rA$p_names)
  rEUniqueList<-rE[-sel,]
  write.table(rEUniqueList, file=paste0("Reactome_unique_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  #sel<-which(rAUniqueList$pvals<0.01) # CAMBIA SOGLIA A PIACIMENTO
  #rAUniqueList[sel,]
  #write.table(rAUniqueList[sel,], file=paste0("Reactome_unique_0.01_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # sel<-which(rEUniqueList$pvals<0.01) # CAMBIA SOGLIA A PIACIMENTO
  # rEUniqueList[sel,]
  # write.table(rEUniqueList[sel,], file=paste0("Reactome_unique_0.01_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # SIGNIFICANCE TEST FOR COMMONS PATHS
  
  common<-which(rA$p_names %in% rE$p_names)
  
  sel<-which(rfA$Var1 %in% rA[common, 'p_names'])
  pNamesCommon<-rfA[sel, 'Var1']
  AFreq<-rfA[sel, 'Freq']
  
  sel<-which(rfE$Var1 %in% rA[common, 'p_names'])
  EFreq<-rfE[sel, 'Freq']
  
  Fmatrix<-data.frame(as.character(pNamesCommon), AFreq, EFreq)
  
  sumA<-sum(Fmatrix$AFreq)
  sumE<-sum(Fmatrix$EFreq)
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATH", "%A", "%E", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    
    posA<-Fmatrix[i, 'AFreq']
    posE<-Fmatrix[i, 'EFreq']
    negA<-sumA-posA
    negE<-sumE-posE
    
    x<-matrix(c(posA, posE, negA, negE), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 1]<-as.character(Fmatrix[i, 1])
      Rmatrix[i, 2]<-posA/sumA*100
      Rmatrix[i, 3]<-posE/sumE*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$QVAL<-p.adjust(Rmatrix$PVAL, method="fdr")
  
  write.table(Rmatrix, file=paste0("Reactome_", index, "_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # GOBP
  
  gobpA<-read.delim(paste0(getwd(), "/GOBP_", index, "_A_0.05.tsv"), header=T, sep="\t", dec=".")
  gobpfA<-read.delim(paste0(getwd(), "/GOBP_", index, "_A_freq.tsv"), header=T, sep="\t", dec=".")
  
  gobpE<-read.delim(paste0(getwd(), "/GOBP_", index, "_E_0.05.tsv"), header=T, sep="\t", dec=".")
  gobpfE<-read.delim(paste0(getwd(), "/GOBP_", index, "_E_freq.tsv"), header=T, sep="\t", dec=".")
  
  sel<-which(gobpA$p_names %in% gobpE$p_names)
  gobpAUniqueList<-gobpA[-sel,]
  sel<-which(gobpE$p_names %in% gobpA$p_names)
  gobpEUniqueList<-gobpE[-sel,]
  
  write.table(gobpAUniqueList, file=paste0("GOBP_unique_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  write.table(gobpEUniqueList, file=paste0("GOBP_unique_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # sel<-which(gobpAUniqueList$pvals<0.05)
  # gobpAUniqueList[sel,]
  # write.table(gobpAUniqueList[sel,], file=paste0("GOBP_unique_0.05_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  # 
  # sel<-which(gobpEUniqueList$pvals<0.05)
  # gobpEUniqueList[sel,]
  # write.table(gobpEUniqueList[sel,], file=paste0("GOBP_unique_0.05_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  
  
  # SIGNIFICANCE TEST FOR COMMONS PATHS
  
  common<-which(gobpA$p_names %in% gobpE$p_names)
  
  sel<-which(gobpfA$Var1 %in% gobpA[common, 'p_names'])
  pNamesCommon<-gobpfA[sel, 'Var1']
  AFreq<-gobpfA[sel, 'Freq']
  
  sel<-which(gobpfE$Var1 %in% gobpA[common, 'p_names'])
  EFreq<-gobpfE[sel, 'Freq']
  
  Fmatrix<-data.frame(as.character(pNamesCommon), AFreq, EFreq)
  
  sumA<-sum(Fmatrix$AFreq)
  sumE<-sum(Fmatrix$EFreq)
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATH", "%A", "%E", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    
    posA<-Fmatrix[i, 'AFreq']
    posE<-Fmatrix[i, 'EFreq']
    negA<-sumA-posA
    negE<-sumE-posE
    
    x<-matrix(c(posA, posE, negA, negE), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 1]<-as.character(Fmatrix[i, 1])
      Rmatrix[i, 2]<-posA/sumA*100
      Rmatrix[i, 3]<-posE/sumE*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$QVAL<-p.adjust(Rmatrix$PVAL, method="fdr")
  
  write.table(Rmatrix, file=paste0("GOBP_", index, "_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  
  # GOMF
  
  gomfA<-read.delim(paste0(getwd(), "/GOMF_", index, "_A_0.05.tsv"), header=T, sep="\t", dec=".")
  gomffA<-read.delim(paste0(getwd(), "/GOMF_", index, "_A_freq.tsv"), header=T, sep="\t", dec=".")
  
  gomfE<-read.delim(paste0(getwd(), "/GOMF_", index, "_E_0.05.tsv"), header=T, sep="\t", dec=".")
  gomffE<-read.delim(paste0(getwd(), "/GOMF_", index, "_E_freq.tsv"), header=T, sep="\t", dec=".")
  
  sel<-which(gomfA$p_names %in% gomfE$p_names)
  gomfAUniqueList<-gomfA[-sel,]
  sel<-which(gomfE$p_names %in% gomfA$p_names)
  gomfEUniqueList<-gomfE[-sel,]
  
  write.table(gomfAUniqueList, file=paste0("GOMF_unique_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  write.table(gomfEUniqueList, file=paste0("GOMF_unique_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # sel<-which(gomfAUniqueList$pvals<0.01)
  # gomfAUniqueList[sel,]
  # write.table(gomfAUniqueList[sel,], file=paste0("GOMF_unique_0.01_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  # 
  # sel<-which(gomfEUniqueList$pvals<0.01)
  # gomfEUniqueList[sel,]
  # write.table(gomfEUniqueList[sel,], file=paste0("GOMF_unique_0.01_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  
  # SIGNIFICANCE TEST FOR COMMONS PATHS
  
  common<-which(gomfA$p_names %in% gomfE$p_names)
  
  sel<-which(gomffA$Var1 %in% gomfA[common, 'p_names'])
  pNamesCommon<-gomffA[sel, 'Var1']
  AFreq<-gomffA[sel, 'Freq']
  
  sel<-which(gomffE$Var1 %in% gomfA[common, 'p_names'])
  EFreq<-gomffE[sel, 'Freq']
  
  Fmatrix<-data.frame(as.character(pNamesCommon), AFreq, EFreq)
  
  sumA<-sum(Fmatrix$AFreq)
  sumE<-sum(Fmatrix$EFreq)
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATH", "%A", "%E", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    
    posA<-Fmatrix[i, 'AFreq']
    posE<-Fmatrix[i, 'EFreq']
    negA<-sumA-posA
    negE<-sumE-posE
    
    x<-matrix(c(posA, posE, negA, negE), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 1]<-as.character(Fmatrix[i, 1])
      Rmatrix[i, 2]<-posA/sumA*100
      Rmatrix[i, 3]<-posE/sumE*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=5)
  
  write.table(Rmatrix, file=paste0("GOMF_", index, "_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # GOCC
  
  goccA<-read.delim(paste0(getwd(), "/GOCC_", index, "_A_0.05.tsv"), header=T, sep="\t", dec=".")
  goccfA<-read.delim(paste0(getwd(), "/GOCC_", index, "_A_freq.tsv"), header=T, sep="\t", dec=".")
  
  goccE<-read.delim(paste0(getwd(), "/GOCC_", index, "_E_0.05.tsv"), header=T, sep="\t", dec=".")
  goccfE<-read.delim(paste0(getwd(), "/GOCC_", index, "_E_freq.tsv"), header=T, sep="\t", dec=".")
  
  sel<-which(goccA$p_names %in% goccE$p_names)
  goccAUniqueList<-goccA[-sel,]
  sel<-which(goccE$p_names %in% goccA$p_names)
  goccEUniqueList<-goccE[-sel,]
  
  write.table(goccAUniqueList, file=paste0("GOCC_unique_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  write.table(goccEUniqueList, file=paste0("GOCC_unique_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # sel<-which(goccAUniqueList$pvals<0.01)
  # goccAUniqueList[sel,]
  # write.table(goccAUniqueList[sel,], file=paste0("GOCC_unique_0.01_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  # 
  # sel<-which(goccEUniqueList$pvals<0.01)
  # goccEUniqueList[sel,]
  # write.table(goccEUniqueList[sel,], file=paste0("GOCC_unique_0.01_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  
  # SIGNIFICANCE TEST FOR COMMONS PATHS
  
  common<-which(goccA$p_names %in% goccE$p_names)
  
  sel<-which(goccfA$Var1 %in% goccA[common, 'p_names'])
  pNamesCommon<-goccfA[sel, 'Var1']
  AFreq<-goccfA[sel, 'Freq']
  
  sel<-which(goccfE$Var1 %in% goccA[common, 'p_names'])
  EFreq<-goccfE[sel, 'Freq']
  
  Fmatrix<-data.frame(as.character(pNamesCommon), AFreq, EFreq)
  
  sumA<-sum(Fmatrix$AFreq)
  sumE<-sum(Fmatrix$EFreq)
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATH", "%A", "%E", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    
    posA<-Fmatrix[i, 'AFreq']
    posE<-Fmatrix[i, 'EFreq']
    negA<-sumA-posA
    negE<-sumE-posE
    
    x<-matrix(c(posA, posE, negA, negE), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 1]<-as.character(Fmatrix[i, 1])
      Rmatrix[i, 2]<-posA/sumA*100
      Rmatrix[i, 3]<-posE/sumE*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=5)
  
  write.table(Rmatrix, file=paste0("GOCC_", index, "_comparison.tsv", sep="\t"), row.names=F, col.names=T, dec=".", quote=F)
  
  
  ## KEGG
  
  keggA<-read.delim(paste0(getwd(), "/KEGG_", index, "_A_0.05.tsv"), header=T, sep="\t", dec=".")
  keggfA<-read.delim(paste0(getwd(), "/KEGG_", index, "_A_freq.tsv"), header=T, sep="\t", dec=".")
  
  keggE<-read.delim(paste0(getwd(), "/KEGG_", index, "_E_0.05.tsv"), header=T, sep="\t", dec=".")
  keggfE<-read.delim(paste0(getwd(), "/KEGG_", index, "_E_freq.tsv"), header=T, sep="\t", dec=".")
  
  sel<-which(keggA$p_names %in% keggE$p_names)
  keggAUniqueList<-keggA[-sel,]
  sel<-which(keggE$p_names %in% keggA$p_names)
  keggEUniqueList<-keggE[-sel,]
  
  write.table(keggAUniqueList, file=paste0("KEGG_unique_", index, "_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  write.table(keggEUniqueList, file=paste0("KEGG_unique_", index, "_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  # sel<-which(keggAUniqueList$pvals<0.01)
  # keggAUniqueList[sel,]
  # write.table(keggAUniqueList[sel,], file=paste0("KEGG_", index, "_unique_0.01_A.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  # 
  # sel<-which(keggEUniqueList$pvals<0.01)
  # keggEUniqueList[sel,]
  # write.table(keggEUniqueList[sel,], file=paste0("KEGG_", index, "_unique_0.01_E.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  
  # SIGNIFICANCE TEST FOR COMMONS PATHS
  
  common<-which(keggA$p_names %in% keggE$p_names)
  
  sel<-which(keggfA$Var1 %in% keggA[common, 'p_names'])
  pNamesCommon<-keggfA[sel, 'Var1']
  AFreq<-keggfA[sel, 'Freq']
  
  sel<-which(keggfE$Var1 %in% keggA[common, 'p_names'])
  EFreq<-keggfE[sel, 'Freq']
  
  Fmatrix<-data.frame(as.character(pNamesCommon), AFreq, EFreq)
  
  sumA<-sum(Fmatrix$AFreq)
  sumE<-sum(Fmatrix$EFreq)
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATH", "%A", "%E", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    
    posA<-Fmatrix[i, 'AFreq']
    posE<-Fmatrix[i, 'EFreq']
    negA<-sumA-posA
    negE<-sumE-posE
    
    x<-matrix(c(posA, posE, negA, negE), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 1]<-as.character(Fmatrix[i, 1])
      Rmatrix[i, 2]<-posA/sumA*100
      Rmatrix[i, 3]<-posE/sumE*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=5)
  
  write.table(Rmatrix, file=paste0("KEGG_", index, "_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
              
}


#---------------------------------------------------------------------------------------------------------------
#-------------------------------------- FIND SIGNIFICANT GENES BETWEEN GROUPS ----------------------------------
#---------------------------------------------------------------------------------------------------------------

find_genes_sig<-function(caseTab, ctrlTab, event){
  geneCase<-read.delim(paste0(getwd(), "/", caseTab), header=T, sep="\t", dec=".")
  geneCtrl<-read.delim(paste0(getwd(), "/", ctrlTab), header=T, sep="\t", dec=".")
  
  mgeneCase<-as.matrix(geneCase[,2:dim(geneCase)[2]])
  mgeneCtrl<-as.matrix(geneCtrl[,2:dim(geneCtrl)[2]])  
  
  sel<-which(mgeneCase>1)
  if(length(sel)>1) { mgeneCase[sel]=1 }
  
  sel<-which(mgeneCtrl>1)
  if(length(sel)>1) { mgeneCtrl[sel]=1 }
  
  geneCase[,2:dim(geneCase)[2]]<-mgeneCase
  geneCtrl[,2:dim(geneCtrl)[2]]<-mgeneCtrl
  
  rm(mgeneCase)
  rm(mgeneCtrl)
  
  caseCounts<-apply(geneCase[,2:dim(geneCase)[2]], 2, sum)
  ctrlCounts<-apply(geneCtrl[,2:dim(geneCtrl)[2]], 2, sum)
  
  Fmatrix<-data.frame(caseCounts, ctrlCounts)
  
  sel<-grep("^LOC", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("^LINC", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("orf", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("^TTTY", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("SNORD", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("^FAM", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("KRT", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("RYR", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("PCLO", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("CSMD", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("CNTN", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  sel<-grep("KIA", rownames(Fmatrix))
  if(length(sel)>0) Fmatrix<-Fmatrix[-sel,]
  
  sumCase<-dim(geneCase)[1]
  sumCtrl<-dim(geneCtrl)[1]
  
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("GENE", "%Case", "%Ctrl", "PVAL", "QVAL")
  
  for( i in 1:dim(Fmatrix)[1]){
    cat(i, "\n")
    posCase<-Fmatrix[i, 'caseCounts']
    posCtrl<-Fmatrix[i, 'ctrlCounts']
    negCase<-sumCase-posCase
    negCtrl<-sumCtrl-posCtrl
    
    x<-matrix(c(posCase, posCtrl, negCase, negCtrl), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
      
      Rmatrix[i, 2]<-posCase/sumCase*100
      Rmatrix[i, 3]<-posCtrl/sumCtrl*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$GENE<-rownames(Fmatrix)
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=5)
  
  write.table(Rmatrix, file=paste0(event,"_genes_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  return(Rmatrix)
}

#m<-find_genes_sig("delCountsPerPatientTableCase.tsv", "delCountsPerPatientTableCtrl.tsv", "del")
#m<-find_genes_sig("dupCountsPerPatientTableCase.tsv", "dupCountsPerPatientTableCtrl.tsv", "dup")
m<-find_genes_sig("gainAllCountsPerPatientTableA.tsv", "gainAllCountsPerPatientTableE.tsv", "gain")
m<-find_genes_sig("lossAllCountsPerPatientTableA.tsv", "lossAllCountsPerPatientTableE.tsv", "loss")

# geneCase<-read.delim(paste0(getwd(), "/lossGenesCase.tsv"), header=T, sep="\t", dec=".")
# geneCtrl<-read.delim(paste0(getwd(), "/lossGenesCtrl.tsv"), header=T, sep="\t", dec=".")
# 
# sel<-grep("^LOC", geneCase$GENE)
# if(length(sel)>0) geneCase<-geneCase[-sel,]
# sel<-grep("^LINC", geneCase$GENE)
# if(length(sel)>0) geneCase<-geneCase[-sel,]
# sel<-grep("orf", geneCase$GENE)
# if(length(sel)>0) geneCase<-geneCase[-sel,]
# 
# sel<-grep("^LOC", geneCtrl$GENE)
# if(length(sel)>0) geneCtrl<-geneCtrl[-sel,]
# sel<-grep("^LINC", geneCtrl$GENE)
# if(length(sel)>0) geneCtrl<-geneCtrl[-sel,]
# sel<-grep("orf", geneCtrl$GENE)
# if(length(sel)>0) geneCtrl<-geneCtrl[-sel,]
# 
# 
# Fmatrix<-data.frame(geneCase$GENE, geneCase$COUNTS, geneCtrl$COUNTS)
# 
# uniqueGenesA<-Fmatrix[which(Fmatrix[,2]!=0 & Fmatrix[,3]==0),]
# uniqueGenesE<-Fmatrix[which(Fmatrix[,2]==0 & Fmatrix[,3]!=0),]
# Fmatrix<-Fmatrix[which(Fmatrix[,2]!=0 & Fmatrix[,3]!=0),]
# 
# sumA<-sum(Fmatrix$geneCase.COUNTS)
# sumE<-sum(Fmatrix$geneCtrl.COUNTS)
# 
# Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
# colnames(Rmatrix)<-c("GENE", "%A", "%E", "PVAL", "QVAL")
# 
# for( i in 1:dim(Fmatrix)[1]){
#   
#   posA<-Fmatrix[i, 'geneCase.COUNTS']
#   posE<-Fmatrix[i, 'geneCtrl.COUNTS']
#   negA<-sumA-posA
#   negE<-sumE-posE
#   
#   x<-matrix(c(posA, posE, negA, negE), ncol=2, byrow=T)
#   if(sum(!is.na(x))==4){
#     
#     pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=5)
#     
#     Rmatrix[i, 1]<-as.character(Fmatrix[i, 1])
#     Rmatrix[i, 2]<-posA/sumA*100
#     Rmatrix[i, 3]<-posE/sumE*100
#     Rmatrix[i, 4]<-pval
#   } 
# }
# Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=5)
# 
# write.table(Rmatrix, file="loss_genes_comparison.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)
# write.table(uniqueGenesA, file="loss_genes_unique_A.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)
# write.table(uniqueGenesE, file="loss_genes_unique_E.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)


#---------------------------------------------------------------------------------------------------------------

# ADDING GENE NAMES TO REACTOME TABLE

reactome_files = list.files(getwd(), "geneAnnoReactome")
for(i in 1:length(reactome_files)){
  file = reactome_files[i]
  reactome = read.delim(file, stringsAsFactors = F)
  symbol = unlist(mapIds(org.Hs.eg.db, keys=as.character(reactome$ENTREZID), column="SYMBOL", keytype="ENTREZID"))
  reactome$GENE = symbol
  write.table(reactome, file=paste0(file), sep="\t", row.names=F, col.names=T, quote=F)
}


# geneAnnoReactome_E.somatic_variants_filtered_rev <- read.delim("/media/marco/DATA2/NGS/WES/Projects/Results_AE/geneAnnoReactome_E.somatic_variants_filtered_rev.tsv")
# which(is.na(geneAnnoReactome_E.somatic_variants_filtered_rev$ENTREZID))
# geneAnnoReactome_E.somatic_variants_filtered_rev<-geneAnnoReactome_E.somatic_variants_filtered_rev[-2460,]
# SYMBOL<-unlist(mapIds(org.Hs.eg.db, keys=as.character(geneAnnoReactome_E.somatic_variants_filtered_rev$ENTREZID), column="SYMBOL", keytype="ENTREZID"))
# length(SYMBOL)==dim(geneAnnoReactome_E.somatic_variants_filtered_rev)[1]
# geneAnnoReactome_E.somatic_variants_filtered_rev<-cbind(geneAnnoReactome_E.somatic_variants_filtered_rev, SYMBOL)
# write.table(geneAnnoReactome_E.somatic_variants_filtered_rev, file="geneAnnoReactome_E.somatic_variants_filtered_rev.tsv", sep="\t", row.names=F, col.names=T, quote=F)
# 
# geneAnnoReactome_A.somatic_variants_filtered_rev <- read.delim("/media/marco/DATA2/NGS/WES/Projects/Results_AE/geneAnnoReactome_A.somatic_variants_filtered_rev.tsv")
# which(is.na(geneAnnoReactome_A.somatic_variants_filtered_rev$ENTREZID))
# #geneAnnoReactome_A.somatic_variants_filtered_rev<-geneAnnoReactome_A.somatic_variants_filtered_rev[-2460,]
# SYMBOL<-unlist(mapIds(org.Hs.eg.db, keys=as.character(geneAnnoReactome_A.somatic_variants_filtered_rev$ENTREZID), column="SYMBOL", keytype="ENTREZID"))
# length(SYMBOL)==dim(geneAnnoReactome_A.somatic_variants_filtered_rev)[1]
# geneAnnoReactome_A.somatic_variants_filtered_rev<-cbind(geneAnnoReactome_A.somatic_variants_filtered_rev, SYMBOL)
# write.table(geneAnnoReactome_A.somatic_variants_filtered_rev, file="geneAnnoReactome_A.somatic_variants_filtered_rev.tsv", sep="\t", row.names=F, col.names=T, quote=F)

#---------------------------------------------------------------------------------------------------------------



