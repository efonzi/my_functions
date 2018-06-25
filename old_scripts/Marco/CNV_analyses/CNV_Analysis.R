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

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(MASS)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTIONS

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createRepTable<-function(repTab){
  
  chr<-matrix(unlist(strsplit(as.character(repTab[,3]), split=c(":"), fixed=T)), ncol=2, byrow=T)
  
  chr[,1]<-gsub("chr", "", chr[,1])
  
  chr[which(chr[, 1]=="X"), 1]<-23
  
  chr[which(chr[, 1]=="Y"), 1]<-24
  
  pos<-matrix(unlist(strsplit(as.character(chr[,2]), split=c("-"), fixed=T)), ncol=2, byrow=T)
  
  pos[,1]<-gsub(",", "", pos[, 1])
  
  pos[,2]<-gsub(",", "", pos[, 2])
  
  cnv<-repTab[,4]
  
  cnv[cnv=="CN Gain"]<-"gain"
  
  cnv[cnv=="CN Loss"]<-"loss"
  
  cnv[cnv=="High Copy Gain"]<-"duplication"
  
  cnv[cnv=="Homozygous Copy Loss"]<-"deletion"
  
  cnv[cnv=="LOH"]<-"UPD"
  
  repTab<-data.frame(as.numeric(chr[,1]), as.numeric(pos[,1]), as.numeric(pos[,2]), repTab$Probe.Median, cnv, repTab$Sample, repTab$Count.of.Gene.Symbols, repTab$Gene.Symbol, stringsAsFactors=F)
  
  colnames(repTab)<-c("CHR", "START", "END", "MEDIAN", "CNV_TYPE", "SAMPLE", "GENE.COUNT", "GENE")
  
  return(repTab)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------




createGeneCountsMatrix<-function(dset){
  
  # Unlist gene names and create a matrix
  
  listed<-dset[which(dset$GENE.COUNT>1),]
  m<-max(dset$GENE.COUNT)
  unlisted<-matrix(rep(NA, m*dim(listed)[1]), nrow=dim(listed)[1], ncol=m)
  for(i in 1:dim(listed)[1]){
    ncol<-listed[i, 'GENE.COUNT']
    #row<-unique(unlist(strsplit(as.character(listed[i, 'GENE']), split=c(","), fixed=T)))
    row<-unlist(lapply(unique(unlist(strsplit(as.character(listed[i, 'GENE']), split=c(","), fixed=T))), function(x) (gsub(" ", "", x))))
    if(length(row)<listed[i, 'GENE.COUNT']){
      ncol<-length(row)
    } 
    tmp<-matrix(row, ncol=ncol, byrow=T)
    unlisted[i, 1:dim(tmp)[2]]<-tmp
  }
  
  # Remove space from genes names
  #unlisted<-apply(unlisted, 2, function(x) (gsub(" ", "", x)))
  
  # Create a matrix with single genes
  genes<-dset[which(dset$GENE.COUNT==1),]
  m<-dim(unlisted)[2]
  genesMatrix<-matrix(rep(NA, m*dim(genes)[1]), nrow=dim(genes)[1], ncol=m)
  genesMatrix[,1]<-genes$GENE
  
  if(dim(genesMatrix)[1]>0){
    # Remove space from genes names
    genesMatrix<-apply(genesMatrix, 2, function(x) (gsub(" ", "", x)))
  }
  
  # Bind sample column
  unlisted<-cbind(unlisted, listed$SAMPLE)
  
  if(dim(genesMatrix)[1]>0){
    genesMatrix<-cbind(genesMatrix, genes$SAMPLE)
  
    # Merge matrices
    genesMatrix<-rbind(genesMatrix, unlisted)
    genesMatrix<-data.frame(genesMatrix, stringsAsFactors=F)
    
  } else {
    
    genesMatrix<-data.frame(unlisted, stringsAsFactors=F)
    
  }
  
  return(genesMatrix)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createNullHyp<-function(dset, col){
  
  R=1000
  n=dim(dset)[1]
  par<-fitdistr(dset[,col], "negative binomial", lower=0.1)
  x<-rnegbin(n, par$estimate[1], par$estimate[2])
  #x<-rpois(n, par$estimate)
  
  nh<-rep(0, n)
  m<-matrix(rep(0, R*n), ncol=R, byrow=F)
    for(i in 1:R){
      
      tmp<-sample(x, n, replace = T)
      
      m[,i]<-tmp
      
    }
  
  nh<-round(apply(m, 1, mean), digits=0)

  return(list(nh=nh, m=m, x=x))
  
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createGeneCountsTable<-function(dset, basename){
  
  geneCountsTable<-createGeneCountsMatrix(dset)
  
  last_col<-dim(geneCountsTable)[2]
  
  geneNamesTmp<-vector()
  geneNames<-vector()
  for(i in 1:(last_col-1)){
    
    geneNamesTmp<-geneCountsTable[!is.na(geneCountsTable[,i]), i]
    
    cat("Genes: ", length(geneNamesTmp), "\t", "Iteration: ", i, "\n")
    
    geneNames<-c(geneNames, geneNamesTmp)
    
    
  }
  
  geneNames<-unique(geneNames)
  
  geneCounts<-rep(0, length(geneNames))
  
  cat("Found ", length(geneNames), " genes.\n")
  
  comp<-function(x, y){
    return(length(which(x==y)))
  }
  
  for(i in 1:length(geneNames)){
    
    cat("Searching gene occurences..\n")
    
    cat("Iteration: ", i, " in ", length(geneNames), "\n")
    
    geneCounts[i]<-comp(geneCountsTable, geneNames[i])
    
  }
  
  geneNames<-sapply(geneNames, function(x){as.character(x)})
  geneCounts<-sapply(geneCounts, function(x){as.numeric(x)})
  
  CountsTable<-data.frame(geneNames, geneCounts)
  colnames(CountsTable)=c("GENE", "COUNTS")
  

  ### Remove genes from geneCountsTable (needed for geneVsPatients table)
  
  cat("Remove genes from geneCountsTable (needed for geneVsPatients table)\n")
  
  for(j in 1:(last_col-1)){
    sel<-which(geneCountsTable[, j] %in% CountsTable[,1])
    if(length(sel)>0) geneCountsTable[-sel, j]="NA"
    else geneCountsTable[,j]="NA"
  }
  
  
  ### Shift left genes
  
  cat("Shift left genes\n")
  
  m<-dim(geneCountsTable)[1]
  n<-dim(geneCountsTable)[2]
  
  tmpMat<-matrix(rep("NA", m*n), ncol=n)
  
  for(i in 1:m){
    k=1
    for(j in 1:(n-1)){
      if(geneCountsTable[i,j]!="NA") {
        tmpMat[i,k]<-geneCountsTable[i,j]
        k=k+1
      }
    }   
  }
  tmpMat[,n]<-geneCountsTable[,n]
  geneCountsTable<-tmpMat
  sel<-which(geneCountsTable[,1]=="NA")
  if(length(sel>0)) geneCountsTable<-geneCountsTable[-sel,]  
  
  ### Resize geneCountsTable
  
  cat("Resize geneCountsTable\n")
  
  m=0
  for(i in 1:(last_col-1)){
    sel<-which(geneCountsTable[, i]!="NA")
    if(length(sel)>0) m=m+1
  }
  geneCountsTable<-data.frame(geneCountsTable[, 1:m], geneCountsTable[, last_col], stringsAsFactors=F)
  sel<-which(geneCountsTable[,1]=="NA")
  if(length(sel)>0) geneCountsTable<-geneCountsTable[-sel,]
  
  ### Gene per patients table
  
  cat("Gene per patients table\n")
  
  last_col<-m+1
  pats<-levels(as.factor(geneCountsTable[, last_col]))
  m<-length(pats)
  
  geneNamesTmp<-vector()
  geneNames<-vector()
  
  for(i in 1:(last_col-1)){
    geneNamesTmp<-geneCountsTable[,i]
    geneNames<-c(geneNames, geneNamesTmp)
  }
  
  geneNames<-unique(geneNames)
  sel<-which(geneNames=="NA")
  if(length(sel)>0) geneNames<-geneNames[-sel]
  n<-length(geneNames)
  
  mat<-matrix(rep(0, m*n), ncol=n, byrow=F)
  
  CountsPerPatientsTable<-data.frame(pats, mat, stringsAsFactors=F)
  colnames(CountsPerPatientsTable)<-c("Sample", geneNames)
  
  for(i in 1:dim(geneCountsTable)[1]){
    for(j in 1:(last_col-1)){
      row=as.character(geneCountsTable[i, last_col])
      rowSel<-which(CountsPerPatientsTable$Sample==row)
      col=as.character(geneCountsTable[i, j])
      if(col!="NA"){
        colSel<-which(colnames(CountsPerPatientsTable)==col)
        count=CountsPerPatientsTable[rowSel, colSel]
        count=count+1
        CountsPerPatientsTable[rowSel, colSel]=count
      }
    }
  }
  
  rm(tmpMat)
  
#   tmpTable<-(table(geneCountsTable[,1], geneCountsTable[,last_col]))
#   last_col<-dim(tmpTable)[2]
#   CountsPerPatientTable<-data.frame(row.names(tmpTable), tmpTable[,1:last_col])
#   names_vector<-colnames(tmpTable)
#   names_vector<-c("GENE", names_vector)
#   names(CountsPerPatientTable)<-names_vector
  
  #write.table(CountsTable, file=paste(basename, "CountsTable.tsv", sep=""), sep="\t", row.names=F, col.names=T)
  write.table(CountsPerPatientsTable, file=paste0("arrayTable_", basename, ".tsv"), sep="\t", row.names=F, col.names=T)
  
  
  return(list(CountsTable=CountsTable, CountsPerPatientsTable=CountsPerPatientsTable))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createGeneCountsTableThrs<-function(dset, alpha, basename){
  
  # Bootstrapped
  
  
  ### Generate Counts table
  
  geneCountsTable<-createGeneCountsMatrix(dset)
  
  
  last_col<-dim(geneCountsTable)[2]
  
  geneNamesTmp<-vector()
  geneNames<-vector()
  for(i in 1:(last_col-1)){
    
    geneNamesTmp<-geneCountsTable[!is.na(geneCountsTable[,i]), i]
    
    cat("Genes: ", length(geneNamesTmp), "\t", "Iteration: ", i, "\n")
    
    geneNames<-c(geneNames, geneNamesTmp)
    
    
  }
  
  geneNames<-unique(geneNames)
  
  geneCounts<-rep(0, length(geneNames))
  
  for(i in 1:length(geneNames)){
    
    cat("Gene: ", geneNames[i], "\t", "Iteration: ", i, "\n")
    
    geneCounts[i]<-length(which(geneCountsTable==geneNames[i]))
    
  }
  
  geneNames<-sapply(geneNames, function(x){as.character(x)})
  geneCounts<-sapply(geneCounts, function(x){as.numeric(x)})
  
  CountsTable<-data.frame(geneNames, geneCounts, stringsAsFactors=F)
  colnames(CountsTable)=c("GENE", "COUNTS")
  
  
  ### Bootstrap
  cat("Bootstrap\n")
  
  nh<-createNullHyp(CountsTable, 2)
  
  m<-dim(CountsTable)[1]
  
  pvals<-rep(0, m)
  
  th0<-sum(sapply(nh[1], function(x){as.numeric(x)}))
  th1<-length(levels(as.factor(CountsTable$GENE)))
  
  for(i in 1:m){
    
    c0<-(sapply(nh[1], function(x){as.numeric(x)}))[i]
    c1<-CountsTable[i, 'COUNTS']
    ct<-matrix(c(c0, th0-c0, c1, th1-c1), nrow=2, byrow=F)
    ft<-fisher.test(ct)
    #csq<-chisq.test(ct, simulate.p.value=T)
    
    pvals[i]<-ft$p.value
    
  }
  
  CountsTable<-data.frame(CountsTable, pvals, stringsAsFactors=F)
  sel<-which(CountsTable$pvals<alpha)
  CountsTable<-CountsTable[sel,]  
  
  ### Remove genes from geneCountsTable (needed for geneVsPatients table)
  
  cat("Remove genes from geneCountsTable (needed for geneVsPatients table)\n")
  
  for(j in 1:(last_col-1)){
    sel<-which(geneCountsTable[, j] %in% CountsTable[,1])
    if(length(sel)>0) geneCountsTable[-sel, j]="NA"
       else geneCountsTable[,j]="NA"
  }
  
  
  ### Shift left genes
  
  cat("Shift left genes\n")
  
  m<-dim(geneCountsTable)[1]
  n<-dim(geneCountsTable)[2]
  
  tmpMat<-matrix(rep("NA", m*n), ncol=n)
  
  for(i in 1:m){
    k=1
    for(j in 1:(n-1)){
      if(geneCountsTable[i,j]!="NA") {
        tmpMat[i,k]<-geneCountsTable[i,j]
        k=k+1
      }
    }   
  }
  tmpMat[,n]<-geneCountsTable[,n]
  geneCountsTable<-tmpMat
  sel<-which(geneCountsTable[,1]=="NA")
  if(length(sel)>0) geneCountsTable<-geneCountsTable[-sel,]  

  ### Resize geneCountsTable
  
  cat("Resize geneCountsTable\n")
  
  m=0
  for(i in 1:(last_col-1)){
    sel<-which(geneCountsTable[, i]!="NA")
    if(length(sel)>0) m=m+1
  }
  geneCountsTable<-data.frame(geneCountsTable[, 1:m], geneCountsTable[, last_col], stringsAsFactors=F)
  sel<-which(geneCountsTable[,1]=="NA")
  if(length(sel)>0) geneCountsTable<-geneCountsTable[-sel,]
  
  ### Gene per patients table
  
  cat("Gene per patients table\n")
  
  last_col<-m+1
  pats<-levels(as.factor(geneCountsTable[, last_col]))
  m<-length(pats)
  
  geneNamesTmp<-vector()
  geneNames<-vector()
  
  for(i in 1:(last_col-1)){
    geneNamesTmp<-geneCountsTable[,i]
    geneNames<-c(geneNames, geneNamesTmp)
  }
  
  geneNames<-unique(geneNames)
  sel<-which(geneNames=="NA")
  if(length(sel)>0) geneNames<-geneNames[-sel]
  n<-length(geneNames)
  
  mat<-matrix(rep(0, m*n), ncol=n, byrow=F)
  
  CountsPerPatientsTable<-data.frame(pats, mat, stringsAsFactors=F)
  colnames(CountsPerPatientsTable)<-c("Sample", geneNames)
  
  for(i in 1:dim(geneCountsTable)[1]){
    for(j in 1:(last_col-1)){
          row=as.character(geneCountsTable[i, last_col])
          rowSel<-which(CountsPerPatientsTable$Sample==row)
          col=as.character(geneCountsTable[i, j])
          if(col!="NA"){
            colSel<-which(colnames(CountsPerPatientsTable)==col)
            count=CountsPerPatientsTable[rowSel, colSel]
            count=count+1
            CountsPerPatientsTable[rowSel, colSel]=count
          }
    }
  }
  
   
  write.table(CountsTable, file=paste(basename, "CountsTableThrs.tsv", sep=""), sep="\t", row.names=F, col.names=T)
  write.table(CountsPerPatientsTable, file=paste(basename, "CountsPerPatientTableThrs.tsv", sep=""), sep="\t", row.names=F, col.names=T)
  
  rm(tmpMat)
  
  return(list(CountsTableThrs=CountsTable, CountsPerPatientsTableThrs=CountsPerPatientsTable))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createEventsCHRTable<-function(){
  
  dupchr<-table(dupTab$CHR)
  delchr<-table(delTab$CHR)
  gainchr<-table(gainTab$CHR)
  losschr<-table(lossTab$CHR)
  updchr<-table(updTab$CHR)
  m<-matrix(rep(NA, 5*24), ncol=5, byrow=F)
  m[as.numeric(names(dupchr)),1]<-dupchr[names(dupchr)]
  m[as.numeric(names(delchr)),2]<-delchr[names(delchr)]
  m[as.numeric(names(gainchr)),3]<-gainchr[names(gainchr)]
  m[as.numeric(names(losschr)),4]<-dupchr[names(losschr)]
  m[as.numeric(names(updchr)),5]<-updchr[names(updchr)]
  m[is.na(m)]<-0
  
  chr<-c(1:24)
  m<-cbind(chr, m)
  
  colnames(m)<-c("CHR", "DUPLICATION", "DELETION", "GAIN", "LOSS", "UPD")
  
  write.table(m, file="EventsCountsPerCHRTable.tsv", sep="\t", row.names=F, col.names=T)
  
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createCoOccurenceMatrix<-function(dset){
  
  
  
  return(geneCoOccurrenceMatrix)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


searchMinimalCommonRegion<-function(dset){
  
  # BP TOLERANCE
     
  s1=50000; s2=50000; e1=50000; e2=50000  
  
  chrs<-sort(unique(dset$CHR))
  pats<-unique(dset$SAMPLE)
  events<-dim(dset)[1]
  
  m=events
  n=length(pats)+5
  
  mat<-matrix(rep("NA", m*n), ncol=n)
  
  selected=rep(0, dim(dset)[1])
  tmp<-dset
  
  rowPos=1
  
  for(i in 1:dim(dset)[1]){
    #cat(i, "\n")
    sel<-which((tmp$CHR==dset[i, 'CHR']) & (tmp$START>=dset[i, 'START']-s1) & (tmp$START<=dset[i, 'START']+s2) & (tmp$END>=dset[i, 'END']-e1) & (tmp$END<=dset[i, 'END']+e2))
      
    if(length(sel)>0){
      
      if(length(which(selected[sel]==0))>0){
      
        cat("CHR: ", dset[i, 'CHR'], " START: ", dset[i, 'START'], " END: ", dset[i, 'END'], "\n")
        #cat("m ", m, "\n")
        #cat("n ", n, "\n")
        pats<-unique(dset[sel, 'SAMPLE'])
        mat[rowPos, 1]<-dset[i, 'CHR']
        mat[rowPos, 2]<-dset[i, 'START']
        mat[rowPos, 3]<-dset[i, 'END']
        mat[rowPos, 4]<-dset[i, 'GENE']
        mat[rowPos, 5]<-length(pats)
        col<-6+length(pats)-1
        #cat("col ",col, "\n" )
        mat[rowPos, 6:col]<-pats
        #cat("row ", rowPos, "\n")
        rowPos=rowPos+1
        selected[sel]=1
      }
    }
  }
  
  sel<-which(mat[,1]=="NA")
  if(length(sel)>0) mat<-mat[1:(sel[1]-1), ]
  
  k=0
  for(i in 1:n){
    sel<-which(mat[, i]!="NA")
    if(length(sel)>0) k=k+1
  }
  
  mat<-data.frame(mat[,1:k], stringsAsFactors=F)
  mat[,5]<-as.numeric(mat[,5])
  
  ### Bootstrap
  
  nh<-createNullHyp(mat, 5)
  
  m<-dim(mat)[1]
  
  pvals<-rep(0, m)
  
  th0<-sum(sapply(nh[1], function(x){as.numeric(x)}))
  th1<-sum(mat[,5])
  
  for(i in 1:m){
    
    c0<-(sapply(nh[1], function(x){as.numeric(x)}))[i]
    c1<-mat[i, 5]
    ct<-matrix(c(c0, th0-c0, c1, th1-c1), nrow=2, byrow=F)
    ft<-fisher.test(ct)
    #csq<-chisq.test(ct, simulate.p.value=T)
    
    pvals[i]<-ft$p.value
    
  }
  
  mat<-cbind(mat[,1:5], pvals, mat[, 6:k])
  
  #sel<-which(mat[,5]>1)
  #mat<-mat[sel,]
  
  colnames(mat)<-c("CHR", "START", "END", "GENE", "EVENT_FREQ", "PVAL", "SAMPLE")
  
  return(mat)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


createEventsPatientsTable<-function(repTab){
  
  tmpSet=NULL
  
  patients<-levels(as.factor(repTab$SAMPLE))
  
  for(i in 1:length(patients)){
    
    sel<-which(repTab$SAMPLE==patients[i])
    tmpSet<-repTab[sel,]
    fileName<-patients[i]
    write.table(tmpSet, file=paste(fileName, "SVs.tsv", sep="."), row.names=F, col.names=T, sep="\t", quote=F, dec=".")
  }
  
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("STARTING SCRIPT..\n")

cat("Loading data..\n")

#setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/2017-02-28_enrichmentdatinonfiltrati")

REPFILENAME="dgv_filtered.tsv"

# LOAD REPORT FILE

#dset<-read.delim(paste(getwd(), REPFILENAME, sep="/") , header=T, sep="\t", stringsAsFactors=F)
dset<-read.delim(paste0(getwd(), "/3.txt"), header=T, sep="\t", stringsAsFactors=F)

### EUGENIO
# I added this part because my dataset did not have "GENE.COUNT" column
dset = unique(dset)
dset$GENE.COUNT = 0
dset = dset[dset$GENE != "",]
for(r in 1:dim(dset)[1]){
  genes = dset$GENE[r]
  dset$GENE.COUNT[r] = length(strsplit(genes, ",")[[1]])}

repTab = dset

# del[[2]][,200]
# for(c in 1:dim(del[[2]])[2]){print(3 %in% del[[2]][,c])}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# CREATE REPORT TABS

cat("Report tabs..\n")

#All events
repTab<-createRepTable(dset)

#Gain
gainTab<-repTab[which(repTab$CNV_TYPE=="gain"),]
#sel<-which(gainTab$END-gainTab$START>1000)
#gainTab<-gainTab[sel,]
write.table(gainTab, file="Gain_All.tsv", sep="\t", row.names=F, col.names=T, quote = F)

#Loss
lossTab<-repTab[which(repTab$CNV_TYPE=="loss"),]
#sel<-which(lossTab$END-lossTab$START>1000)
#lossTab<-lossTab[sel,]
write.table(lossTab, file="Loss_All.tsv", sep="\t", row.names=F, col.names=T, quote = F)

#Duplication
dupTab<-repTab[which(repTab$CNV_TYPE=="duplication"),]
write.table(dupTab, file="Dup_All.tsv", sep="\t", row.names=F, col.names=T, quote = F)

#Deletions
delTab<-repTab[which(repTab$CNV_TYPE=="deletion"),]
#sel<-which(delTab$END-delTab$START>1000)
#delTab<-delTab[sel,]
write.table(delTab, file="Del_All.tsv", sep="\t", row.names=F, col.names=T, quote = F)

# #UPD(LOH)
# updTab<-repTab[which(repTab$CNV_TYPE=="UPD"),]
# #sel<-which(updTab$END-updTab$START>1000)
# #updTab<-updTab[sel,]
# write.table(updTab, file="UPD_All.tsv", sep="\t", row.names=F, col.names=T, quote = F)



cat("Counts table..\n")

# GENERATE COUNTS TABLES

# Deletion

del<-createGeneCountsTable(delTab, "del")
delThrs<-createGeneCountsTableThrs(delTab, alpha=0.05, "del")

# Duplication

dup<-createGeneCountsTable(dupTab, "dup")
dupThrs<-createGeneCountsTableThrs(dupTab, alpha=0.05, "dup")

# Gain

gain<-createGeneCountsTable(gainTab, "gain")
gainThrs<-createGeneCountsTableThrs(gainTab, alpha=0.05, "gain")

# Loss

loss<-createGeneCountsTable(lossTab, "loss")
lossThrs<-createGeneCountsTableThrs(lossTab, alpha=0.05, "loss")

# UPD

#upd<-createGeneCountsTable(updTab, "upd")
#updThrs<-createGeneCountsTableThrs(updTab, alpha=0.05, "upd")




#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TABLE EVENTS PER CHR

createEventsCHRTable()

# TABLE EVENTS PER PATIENTS

createEventsPatientsTable(repTab)

# CIRCULAR PLOTS PER PATIENTS



# ANALYZE COMMON REGIONS

delCommon<-searchMinimalCommonRegion(delTab)
write.table(delCommon, file="delMinimalCommonRegion_6.tsv", sep="\t", row.names=F, col.names=T)

dupCommon<-searchMinimalCommonRegion(dupTab)
write.table(dupCommon, file="dupMinimalCommonRegion_6.tsv", sep="\t", row.names=F, col.names=T)

lossCommon<-searchMinimalCommonRegion(lossTab)
write.table(lossCommon, file="lossMinimalCommonRegion_6.tsv", sep="\t", row.names=F, col.names=T)

gainCommon<-searchMinimalCommonRegion(gainTab)
write.table(gainCommon, file="gainMinimalCommonRegion_6.tsv", sep="\t", row.names=F, col.names=T)

#updCommon<-searchMinimalCommonRegion(updTab)
#write.table(updCommon, file="updMinimalCommonRegion_6.tsv", sep="\t", row.names=F, col.names=T)


