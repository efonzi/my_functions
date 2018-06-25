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

# Library

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTIONS

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# getGeneListFromTable<-function(tab){
#   
#   geneList<-unlist(strsplit(tab$genes, ","))
#   dup<-duplicated(geneList)
#   geneList<-geneList[!dup]
#   
#   return(geneList)
#   
# }

filter_table<-function(table){
  sel<-grep("^LOC", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("^LINC", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("orf", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("^TTTY", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("SNORD", colnames(table))
  if(length(sel)>0) table<-table[,-sel]
  sel<-grep("^FAM", colnames(table))
  if(length(sel)>0) table<-table[,-sel]
  sel<-grep("KRT", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("RYR", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("PCLO", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("CSMD", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("CNTN", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  sel<-grep("KIA", colnames(table))
  if(length(sel)>0) table<-table[, -sel]
  return(table)
}

makelist<-function(tab){
  for(i in 1:dim(selA)[1]){
    k<-which(tab$Sample %in% selA$CHIPID[i])
    if(length(k)>0){
      patA[i]<-tab[k, 'Sample']
    }
  }
  
  for(i in 1:dim(selE)[1]){
    k<-which(tab$Sample %in% selE$CHIPID[i])
    if(length(k)>0){
      patE[i]<-tab[k, 'Sample']
    }
  }
  return(list(patA=patA, patE=patE))
}

filterOutGenesChromothripsis <- function(dset){
  
  sel = grep("LOC", colnames(dset))
  sel = sel[sel!=grep("CLOCK", colnames(dset))]
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("LINC", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("^OR[0-9]", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("MUC", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("KRT", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("RYR", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("CSMD", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("NRXN", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("CNTN", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("[.]AS", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("MIR", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("ANK", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("SNORD", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("[.]IT", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  sel = grep("SCARNA", colnames(dset))
  if(length(sel)>0) dset<-dset[, -sel]
  
  return(dset)
}

# Eugenio functions
makelist_wes<-function(tab){
  for(i in 1:dim(selA)[1]){
    k<-which(tab$Sample %in% selA$Sample[i])
    if(length(k)>0){
      patA[i]<-tab[k, 'Sample']
    }
  }
  
  for(i in 1:dim(selE)[1]){
    k<-which(tab$Sample %in% selE$Sample[i])
    if(length(k)>0){
      patE[i]<-tab[k, 'Sample']
    }
  }
  return(list(patA=patA, patE=patE))
}

filter_table_wes<-function(table){
  sel<-grep("^LOC", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("^LINC", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("orf", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("^TTTY", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("SNORD", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("^FAM", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("KRT", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("RYR", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("PCLO", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("CSMD", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("CNTN", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  sel<-grep("KIA", table$GENE)
  if(length(sel)>0) table<-table[-sel,]
  return(table)
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

genesCaseCtrlComparison <- function(geneTable, selCase, selCtrl){
  
  geneCase = geneTable[which(geneTable$Sample %in% selCase$Sample),]
  geneCtrl = geneTable[which(geneTable$Sample %in% selCtrl$Sample),]
  
  caseCounts<-apply(geneCase[,2:dim(geneCase)[2]], 2, sum)
  ctrlCounts<-apply(geneCtrl[,2:dim(geneCtrl)[2]], 2, sum)
  
  Fmatrix<-data.frame(caseCounts, ctrlCounts)
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
  return(Rmatrix)
  
}


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# CURRENT_DIR="/media/DATA1/DATA/SNP/RAW/cnv_analysis/last"
# DGV_DIR="/media/DATA1/DATA/Genomes/DGV"
# 
# setwd(CURRENT_DIR)

# DGV files are located in /media/DATA1/DATA/Genomes/DGV

# delDGV_GRCh37<-read.delim("/media/DATA1/DATA/Genomes/DGV/delDGV_GRCh37.tsv", stringsAsFactors=F)
# genes<-getGeneListFromTable(delDGV_GRCh37)
# write.table(genes, file="del_DGV_genes_list.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# 
# dupDGV_GRCh37<-read.delim("/media/DATA1/DATA/Genomes/DGV/dupDGV_GRCh37.tsv", stringsAsFactors=F)
#genes<-getGeneListFromTable(dupDGV_GRCh37)
#write.table(genes, file="dup_DGV_genes_list.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

#gainDGV_GRCh37<-read.delim("/media/DATA1/DATA/Genomes/DGV/gainDGV_GRCh37.tsv", stringsAsFactors=F)
#genes<-getGeneListFromTable(gainDGV_GRCh37)
#write.table(genes, file="gain_DGV_genes_list.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

#lossDGV_GRCh37<-read.delim("/media/DATA1/DATA/Genomes/DGV/lossDGV_GRCh37.tsv", stringsAsFactors=F)
#genes<-getGeneListFromTable(lossDGV_GRCh37)
#write.table(genes, file="loss_DGV_genes_list.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# delList<-read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/del_DGV_genes_list_updated_cleaned.tsv", stringsAsFactors=FALSE)
# delCountsTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/delCountsTable.tsv", stringsAsFactors=FALSE)
# sel<-which(delCountsTable$GENE %in% delList[,1])
# removed<-delCountsTable[sel,]
# filtered<-delCountsTable[-sel,]
# write.table(removed, file="del_DGV_removed_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(filtered, file="del_DGV_filtered_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# 
# dupList<-read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/dup_DGV_genes_list_updated_cleaned.tsv", stringsAsFactors=FALSE)
# dupCountsTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/dupCountsTable.tsv", stringsAsFactors=FALSE)
# sel<-which(dupCountsTable$GENE %in% dupList[,1])
# removed<-dupCountsTable[sel,]
# filtered<-dupCountsTable[-sel,]
# write.table(removed, file="dup_DGV_removed_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(filtered, file="dup_DGV_filtered_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# 
# lossList<-read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/loss_DGV_genes_list_updated_cleaned.tsv", stringsAsFactors=FALSE)
# lossCountsTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/lossCountsTable.tsv", stringsAsFactors=FALSE)
# sel<-which(lossCountsTable$GENE %in% lossList[,1])
# removed<-lossCountsTable[sel,]
# filtered<-lossCountsTable[-sel,]
# write.table(removed, file="loss_DGV_removed_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(filtered, file="loss_DGV_filtered_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# 
# gainList<-read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/gain_DGV_genes_list_updated_cleaned.tsv", stringsAsFactors=FALSE)
# gainCountsTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/gainCountsTable.tsv", stringsAsFactors=FALSE)
# sel<-which(gainCountsTable$GENE %in% gainList[,1])
# removed<-gainCountsTable[sel,]
# filtered<-gainCountsTable[-sel,]
# write.table(removed, file="gain_DGV_removed_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(filtered, file="gain_DGV_filtered_genes.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# 
# ### FILTERING PATIENTS_GENE TABLES
# 
# setwd("/home/mmanfrini/data1/data/Bologna/ProjectAE/last")
# 
# del_DGV_removed_genes <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/del_DGV_removed_genes.tsv", stringsAsFactors=FALSE)
# delCountsPerPatientTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/delCountsPerPatientTable.tsv", stringsAsFactors=FALSE)
# sel<-which(names(delCountsPerPatientTable) %in% del_DGV_removed_genes[,1])
# delCountsPerPatientTable<-delCountsPerPatientTable[, -sel]
# 
# dup_DGV_removed_genes <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/dup_DGV_removed_genes.tsv", stringsAsFactors=FALSE)
# dupCountsPerPatientTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/dupCountsPerPatientTable.tsv", stringsAsFactors=FALSE)
# sel<-which(names(dupCountsPerPatientTable) %in% dup_DGV_removed_genes[,1])
# dupCountsPerPatientTable<-dupCountsPerPatientTable[, -sel]
# 
# gain_DGV_removed_genes <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/gain_DGV_removed_genes.tsv", stringsAsFactors=FALSE)
# gainCountsPerPatientTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/gainCountsPerPatientTable.tsv", stringsAsFactors=FALSE)
# sel<-which(names(gainCountsPerPatientTable) %in% gain_DGV_removed_genes[,1])
# gainCountsPerPatientTable<-gainCountsPerPatientTable[, -sel]
# 
# loss_DGV_removed_genes <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/loss_DGV_removed_genes.tsv", stringsAsFactors=FALSE)
# lossCountsPerPatientTable <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/lossCountsPerPatientTable.tsv", stringsAsFactors=FALSE)
# sel<-which(names(lossCountsPerPatientTable) %in% loss_DGV_removed_genes[,1])
# lossCountsPerPatientTable<-lossCountsPerPatientTable[, -sel]
# 
# write.table(delCountsPerPatientTable, file="delCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(dupCountsPerPatientTable, file="dupCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(gainCountsPerPatientTable, file="gainCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(lossCountsPerPatientTable, file="lossCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

### SNP RECOVERY

# recovery_list<-read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/snp-recovery-list.tsv", stringsAsFactors=FALSE)
# 
# sel<-which(recovery_list$CNV_TYPE=="gain")
# rgains<-recovery_list[sel,]
# 
# sel<-which(recovery_list$CNV_TYPE=="loss")
# rloss<-recovery_list[sel,]
# 
# sel<-which(recovery_list$CNV_TYPE=="duplication")
# rloss<-recovery_list[sel,]
# 
# sel<-which(recovery_list$CNV_TYPE=="deletion")
# rloss<-recovery_list[sel,]
# 
# recover<-function(set, table){
#   for(i in 1:dim(set)[1])
#     if(set[i, 'GENE.COUNT']==1){
#       rkey<-set[i, 'SAMPLE']  #SELEZIONA LA RIGA
#       rsel<-which(table$Sample==rkey) 
#       if(length(rsel)>0) { #LA RIGA ESISTE
#         ckey<-set[i, 'GENE']  #SELEZIONA LA COLONNA
#         csel<-which(colnames(table)==ckey)  
#         if(length(csel)>0) { #LA COLONNA ESISTE
#           table[rsel, csel]<-table[rsel, csel]+1
#         } else { # LA COLONNA NON ESISTE
#           c<-rep(0, dim(table)[1])
#           c[rsel]<-1
#           table<-cbind(table, c)
#           colnames(table)[dim(table)[2]]<-ckey
#         }
#       } else {  #LA RIGA NON ESISTE
#         r<-rep(0, dim(table)[2])
#         table<-rbind(table, r)
#         ckey<-set[i, 'GENE']  #SELEZIONA LA COLONNA
#         csel<-which(colnames(table)==ckey)  
#         if(length(csel)>0) {  #LA COLONNA ESISTE
#           table[dim(table)[1], csel]<-1
#           table[dim(table)[1], 1]<-set[i, 'SAMPLE']
#         } else {  # LA COLONNA NON ESISTE
#           c<-rep(0, dim(table)[1])
#           c[dim(table)[1]]<-1
#           table<-cbind(table, c)
#           table[dim(table)[1], 1]<-set[i, 'SAMPLE']
#           colnames(table)[dim(table)[2]]<-ckey
#         }
#       }
#     }
#   return(table)
# }
# 
# gainCountsPerPatientTable<-recover(rgains, gainCountsPerPatientTable)
# lossCountsPerPatientTable<-recover(rloss, lossCountsPerPatientTable)
# 
# write.table(gainCountsPerPatientTable, file="gainCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
# write.table(lossCountsPerPatientTable, file="lossCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

#----------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------- FOR WES ----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/wes")
wesTable <- unique(read.delim("snv_indel_NOfusion.tsv", stringsAsFactors=FALSE)[,1:9])

# FILTER OUT UNWANTED GENES

wesTable = filter_table_wes(wesTable)

# SEPARATE CASES AND CONTROLS AND MAKE LISTS OF GENES WITH GENE COUNTS
selA<-unique(wesid[which(wesid$A.E=="A"), c(1,4)])
selE<-unique(wesid[which(wesid$A.E=="E"), c(1,4)])

patA<-rep("", dim(selA)[1])
patE<-rep("", dim(selE)[1])

patList<-makelist_wes(wesTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

wesTableA<-wesTable[which(wesTable$Sample %in% patA),]
wesTableE<-wesTable[which(wesTable$Sample %in% patE),]

write.table(wesTableA, file="wesTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(wesTableE, file="wesTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA = as.matrix(table(wesTableA$GENE))
genesE = as.matrix(table(wesTableE$GENE))

write.table(genesA, file="wesGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="wesGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

# TRANSFORM WES VARIANTS TABLE TO PATIENTxGENE MATRIX

patients = sort(unique(wesTable$Sample))
genes = sort(unique(wesTable$GENE))
mat = matrix(rep(0, length(patients)*length(genes)), ncol=(length(genes)))
colnames(mat) = genes
df = data.frame(Sample = patients, mat, stringsAsFactors = F) # aggiungere la colonna "Sample" non sfasa tutto dopo???

for(i in 1:dim(df)[2]){
  gen = colnames(df)[i]
  subset = wesTable[wesTable$GENE == gen,]
  for(j in 1:dim(df)[1]){
    pat = df$Sample[j]
    if(pat %in% subset$Sample){
      pat_tab = table(subset$Sample)
      df[j,i] = as.numeric(pat_tab[pat])
    }
  }
}

#----------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- FOR CNV - 4 tables ----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

### EXTRACTING PER GROUP GENE TABLE

#setwd("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/tables_revised")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/array_AvsE")

#chipid <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/chipnames", stringsAsFactors=FALSE)
chipid <- read.delim("chipnames", stringsAsFactors=FALSE)

selA<-chipid[which(chipid$AvsE=="1" & chipid$CHIPID!=""),]
selE<-chipid[which(chipid$AvsE=="2" & chipid$CHIPID!=""),]

delCountsPerPatientTable <- read.delim("delCountsPerPatientTable_filtered.tsv", sep="\t", stringsAsFactors=FALSE)
delCountsPerPatientTable$Sample<-unlist(lapply(strsplit(delCountsPerPatientTable$Sample, " "), "[[", 1))
dupCountsPerPatientTable <- read.delim("dupCountsPerPatientTable_filtered.tsv", stringsAsFactors=FALSE)
dupCountsPerPatientTable$Sample<-unlist(lapply(strsplit(dupCountsPerPatientTable$Sample, " "), "[[", 1))
gainCountsPerPatientTable <- read.delim("gainCountsPerPatientTable_filtered.tsv", stringsAsFactors=FALSE)
gainCountsPerPatientTable$Sample<-unlist(lapply(strsplit(gainCountsPerPatientTable$Sample, " "), "[[", 1))
lossCountsPerPatientTable <- read.delim("lossCountsPerPatientTable_filtered.tsv", stringsAsFactors=FALSE)
lossCountsPerPatientTable$Sample<-unlist(lapply(strsplit(lossCountsPerPatientTable$Sample, " "), "[[", 1))

### FILTER TABLES FOR NON MEANINGFUL GENES

dupCountsPerPatientTable<-filter_table(dupCountsPerPatientTable)
delCountsPerPatientTable<-filter_table(delCountsPerPatientTable)
gainCountsPerPatientTable<-filter_table(gainCountsPerPatientTable)
lossCountsPerPatientTable<-filter_table(lossCountsPerPatientTable)

### SEPARATE CASES AND CONTROLS, MAKE LISTS OF GENES

patA<-rep("", dim(selA)[1])
patE<-rep("", dim(selE)[1])


### del

patList<-makelist(delCountsPerPatientTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

delCountsPerPatientTableA<-delCountsPerPatientTable[which(delCountsPerPatientTable$Sample %in% patA),]
delCountsPerPatientTableE<-delCountsPerPatientTable[which(delCountsPerPatientTable$Sample %in% patE),]

write.table(delCountsPerPatientTableA, file="delCountsPerPatientTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(delCountsPerPatientTableE, file="delCountsPerPatientTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA<-as.matrix(apply(delCountsPerPatientTableA[, 2:dim(delCountsPerPatientTableA)[2]], 2, sum))
genesE<-as.matrix(apply(delCountsPerPatientTableE[, 2:dim(delCountsPerPatientTableE)[2]], 2, sum))

write.table(genesA, file="delGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="delGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### dup

patList<-makelist(dupCountsPerPatientTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

dupCountsPerPatientTableA<-dupCountsPerPatientTable[which(dupCountsPerPatientTable$Sample %in% patA),]
dupCountsPerPatientTableE<-dupCountsPerPatientTable[which(dupCountsPerPatientTable$Sample %in% patE),]

write.table(dupCountsPerPatientTableA, file="dupCountsPerPatientTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(dupCountsPerPatientTableE, file="dupCountsPerPatientTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA<-as.matrix(apply(dupCountsPerPatientTableA[, 2:dim(dupCountsPerPatientTableA)[2]], 2, sum))
genesE<-as.matrix(apply(dupCountsPerPatientTableE[, 2:dim(dupCountsPerPatientTableE)[2]], 2, sum))

write.table(genesA, file="dupGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="dupGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### gain

patList<-makelist(gainCountsPerPatientTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

gainCountsPerPatientTableA<-gainCountsPerPatientTable[which(gainCountsPerPatientTable$Sample %in% patA),]
gainCountsPerPatientTableE<-gainCountsPerPatientTable[which(gainCountsPerPatientTable$Sample %in% patE),]

write.table(gainCountsPerPatientTableA, file="gainCountsPerPatientTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(gainCountsPerPatientTableE, file="gainCountsPerPatientTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA<-as.matrix(apply(gainCountsPerPatientTableA[, 2:dim(gainCountsPerPatientTableA)[2]], 2, sum))
genesE<-as.matrix(apply(gainCountsPerPatientTableE[, 2:dim(gainCountsPerPatientTableE)[2]], 2, sum))

write.table(genesA, file="gainGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="gainGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### loss

patList<-makelist(lossCountsPerPatientTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

lossCountsPerPatientTableA<-lossCountsPerPatientTable[which(lossCountsPerPatientTable$Sample %in% patA),]
lossCountsPerPatientTableE<-lossCountsPerPatientTable[which(lossCountsPerPatientTable$Sample %in% patE),]

write.table(lossCountsPerPatientTableA, file="lossCountsPerPatientTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(lossCountsPerPatientTableE, file="lossCountsPerPatientTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA<-as.matrix(apply(lossCountsPerPatientTableA[, 2:dim(lossCountsPerPatientTableA)[2]], 2, sum))
genesE<-as.matrix(apply(lossCountsPerPatientTableE[, 2:dim(lossCountsPerPatientTableE)[2]], 2, sum))

write.table(genesA, file="lossGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="lossGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)


#----------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- FOR CNV - 2 tables ----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/array_AvsE")


### MERGE DUP-GAIN AND DEL-LOSS

# DEL-LOSS

del<-read.delim(paste0(getwd(), "/", "delCountsPerPatientTable_filtered.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F)
loss<-read.delim(paste0(getwd(), "/", "lossCountsPerPatientTable_filtered.tsv"), header=T, sep="\t", dec=".",stringsAsFactors = F)

# equals the number of row

sel<-which(loss$Sample %in% del$Sample)
toadd<-loss[-sel, 'Sample']
mtoadd<-data.frame(toadd, matrix(rep(0, length(toadd)*(dim(del)[2]-1)), ncol=(dim(del)[2]-1)))
colnames(mtoadd)<-colnames(del)
del<-rbind(del,mtoadd)

# merge tables

for(i in 2:dim(del)[2]){
  
  sel<-colnames(del)[i]
  if(sel %in% colnames(loss)){
    k<-which(colnames(loss)==sel)
    loss[,k]<-loss[,k]+del[,i]
  }
  else {
    loss<-cbind(loss, del[,i])
    names(loss)[dim(loss)][2]<-names(del)[i]
  }
}

# DUP-GAINS

dup<-read.delim(paste0(getwd(), "/", "dupCountsPerPatientTable_filtered.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F)
gain<-read.delim(paste0(getwd(), "/", "gainCountsPerPatientTable_filtered.tsv"), header=T, sep="\t", dec=".",stringsAsFactors = F)

# equals the numebr of row

sel<-which(gain$Sample %in% dup$Sample)
toadd<-gain[-sel, 'Sample']
mtoadd<-data.frame(toadd, matrix(rep(0, length(toadd)*(dim(dup)[2]-1)), ncol=(dim(dup)[2]-1)))
colnames(mtoadd)<-colnames(dup)
dup<-rbind(dup,mtoadd)

# merge tables

for(i in 2:dim(dup)[2]){
  
  sel<-colnames(dup)[i]
  if(sel %in% colnames(gain)){
    k<-which(colnames(gain)==sel)
    gain[,k]<-gain[,k]+dup[,i]
  }
  else {
    gain<-cbind(gain, dup[,i])
    names(gain)[dim(gain)][2]<-names(dup)[i]
  }
}

write.table(gain, file="gainAllCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(loss, file="lossAllCountsPerPatientTable_filtered.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

### EXTRACTING PER GROUP GENE TABLE

#chipid <- read.delim("/home/mmanfrini/data1/data/Bologna/ProjectAE/last/chipnames", stringsAsFactors=FALSE)
chipid <- read.delim("chipnames", stringsAsFactors=FALSE)

selA<-chipid[which(chipid$AvsE=="1" & chipid$CHIPID!=""),]
selE<-chipid[which(chipid$AvsE=="2" & chipid$CHIPID!=""),]

selA
selE

patA<-rep("", dim(selA)[1])
patE<-rep("", dim(selE)[1])

### gain

gainAllCountsPerPatientTable<-read.delim(paste0(getwd(), "/", "gainAllCountsPerPatientTable_filtered.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F)

patList<-makelist(gainAllCountsPerPatientTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

gainAllCountsPerPatientTableA<-gainAllCountsPerPatientTable[which(gainAllCountsPerPatientTable$Sample %in% patA),]
gainAllCountsPerPatientTableE<-gainAllCountsPerPatientTable[which(gainAllCountsPerPatientTable$Sample %in% patE),]

write.table(gainAllCountsPerPatientTableA, file="gainAllCountsPerPatientTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(gainAllCountsPerPatientTableE, file="gainAllCountsPerPatientTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA<-as.matrix(apply(gainAllCountsPerPatientTableA[, 2:dim(gainAllCountsPerPatientTableA)[2]], 2, sum))
genesE<-as.matrix(apply(gainAllCountsPerPatientTableE[, 2:dim(gainAllCountsPerPatientTableE)[2]], 2, sum))

write.table(genesA, file="gainAllGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="gainAllGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### loss

lossAllCountsPerPatientTable<-read.delim(paste0(getwd(), "/", "lossAllCountsPerPatientTable_filtered.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F)

patList<-makelist(lossAllCountsPerPatientTable)

patA<-unlist(patList[1])
patE<-unlist(patList[2])

patA<-patA[which(patA!="")]
patE<-patE[which(patE!="")]

lossAllCountsPerPatientTableA<-lossAllCountsPerPatientTable[which(lossAllCountsPerPatientTable$Sample %in% patA),]
lossAllCountsPerPatientTableE<-lossAllCountsPerPatientTable[which(lossAllCountsPerPatientTable$Sample %in% patE),]

write.table(lossAllCountsPerPatientTableA, file="lossAllCountsPerPatientTableA.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(lossAllCountsPerPatientTableE, file="lossAllCountsPerPatientTableE.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesA<-as.matrix(apply(lossAllCountsPerPatientTableA[, 2:dim(lossAllCountsPerPatientTableA)[2]], 2, sum))
genesE<-as.matrix(apply(lossAllCountsPerPatientTableE[, 2:dim(lossAllCountsPerPatientTableE)[2]], 2, sum))

write.table(genesA, file="lossAllGenesA.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesE, file="lossAllGenesE.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)


#----------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------- CHROMOTHRIPSIS ----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

#setwd("/home/marco/Dati/Bologna/Bologna/Chromotripsy/last")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis")

dup<-read.delim(paste0(getwd(), "/", "dupCountsPerPatientTable.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F)
gain<-read.delim(paste0(getwd(), "/", "gainCountsPerPatientTable.tsv"), header=T, sep="\t", dec=".",stringsAsFactors = F)
loss<-read.delim(paste0(getwd(), "/", "lossCountsPerPatientTable.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F)
del<-read.delim(paste0(getwd(), "/", "delCountsPerPatientTable.tsv"), header=T, sep="\t", dec=".",stringsAsFactors = F)


## e poi nelle matrici ci sono ancora numeri maggiori di 1!! TRASFORMALE!!

chipid <- read.delim("chipnames.tsv", header=T, sep="\t", stringsAsFactors=FALSE)
chipid$Sample<-unlist(lapply(strsplit(chipid$Sample, " "), "[[", 1))

selCase<-chipid[which(chipid$chromothripsis=="yes" & chipid$Sample!=""),]
selCtrl<-chipid[which(chipid$chromothripsis=="no" & chipid$Sample!=""),]

tablesFiles = list.files(getwd(), "CountsPerPatientTable.tsv")
for(i in 1:length(tablesFiles)){
  
  tableFileName = tablesFiles[i]
  index = strsplit(tableFileName, "CountsPer")[[1]][1]
  
  # import table, filter out unwanted genes, remove "copy" from samples' names, substitute "." with "-" in gene names
  countTable <- read.delim(tableFileName, header=T, sep="\t", stringsAsFactors=FALSE)
  countTable = filterOutGenesChromothripsis(countTable)
  
  countTable$Sample<-unlist(lapply(strsplit(countTable$Sample, " "), "[[", 1))
  colnames(countTable)<-gsub("[.]", "-", colnames(countTable))
  
  genesComparison = genesCaseCtrlComparison(countTable, selCase, selCtrl)
  
  write.table(genesComparison, file=paste0(index, "_genes_comparison.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
}

####################
gain <- read.delim("gainCountsPerPatientTable.tsv", header=T, sep="\t", stringsAsFactors=FALSE)
loss <- read.delim("lossCountsPerPatientTable.tsv", header=T, sep="\t", stringsAsFactors=FALSE)
del <- read.delim("delCountsPerPatientTable.tsv", header=T, sep="\t", stringsAsFactors=FALSE)
dup <- read.delim("dupCountsPerPatientTable.tsv", header=T, sep="\t", stringsAsFactors=FALSE)

table(sam %in% chipid$Sample)
chipid$Sample[which(chipid$Sample %in% sam == FALSE)]
sam = unique(dset$Sample)
sam[which(sam %in% chipid$Sample == FALSE)]
"20090506_172_SNP6.0_(GenomeWideSNP_6)" %in% gain$Sample 
######################

makelist<-function(tab){
  
  k<-which(tab$Sample %in% selCase$SAMPLE)
  if(length(k)>0){
    patCase<-tab[k, 'Sample']
  }
  
  
  k<-which(tab$Sample %in% selCtrl$SAMPLE)
  if(length(k)>0){
    patCtrl<-tab[k, 'Sample']
  }
  
  return(list(patCase=patCase, patCtrl=patCtrl))
}

### del

patCase<-rep("", dim(selCase)[1])
patCtrl<-rep("", dim(selCtrl)[1])

patList<-makelist(delCountsPerPatientTable)

patCase<-unlist(patList[1])
patCtrl<-unlist(patList[2])

patCase<-patCase[which(patCase!="")]
patCtrl<-patCtrl[which(patCtrl!="")]

delCountsPerPatientTableCase<-delCountsPerPatientTable[which(delCountsPerPatientTable$Sample %in% patCase),]
delCountsPerPatientTableCtrl<-delCountsPerPatientTable[which(delCountsPerPatientTable$Sample %in% patCtrl),]

write.table(delCountsPerPatientTableCase, file="delCountsPerPatientTableCase.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(delCountsPerPatientTableCtrl, file="delCountsPerPatientTableCtrl.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesCase<-as.matrix(apply(delCountsPerPatientTableCase[, 2:dim(delCountsPerPatientTableCase)[2]], 2, sum))
genesCtrl<-as.matrix(apply(delCountsPerPatientTableCtrl[, 2:dim(delCountsPerPatientTableCtrl)[2]], 2, sum))

write.table(genesCase, file="delGenesCase.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesCtrl, file="delGenesCtrl.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### dup

patCase<-rep("", dim(selCase)[1])
patCtrl<-rep("", dim(selCtrl)[1])

patList<-makelist(dupCountsPerPatientTable)

patCase<-unlist(patList[1])
patCtrl<-unlist(patList[2])

patCase<-patCase[which(patCase!="")]
patCtrl<-patCtrl[which(patCtrl!="")]

dupCountsPerPatientTableCase<-dupCountsPerPatientTable[which(dupCountsPerPatientTable$Sample %in% patCase),]
dupCountsPerPatientTableCtrl<-dupCountsPerPatientTable[which(dupCountsPerPatientTable$Sample %in% patCtrl),]

write.table(dupCountsPerPatientTableCase, file="dupCountsPerPatientTableCase.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(dupCountsPerPatientTableCtrl, file="dupCountsPerPatientTableCtrl.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesCase<-as.matrix(apply(dupCountsPerPatientTableCase[, 2:dim(dupCountsPerPatientTableCase)[2]], 2, sum))
genesCtrl<-as.matrix(apply(dupCountsPerPatientTableCtrl[, 2:dim(dupCountsPerPatientTableCtrl)[2]], 2, sum))

write.table(genesCase, file="dupGenesCase.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesCtrl, file="dupGenesCtrl.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### gain

patCase<-rep("", dim(selCase)[1])
patCtrl<-rep("", dim(selCtrl)[1])

patList<-makelist(gainCountsPerPatientTable)

patCase<-unlist(patList[1])
patCtrl<-unlist(patList[2])

patCase<-patCase[which(patCase!="")]
patCtrl<-patCtrl[which(patCtrl!="")]

gainCountsPerPatientTableCase<-gainCountsPerPatientTable[which(gainCountsPerPatientTable$Sample %in% patCase),]
gainCountsPerPatientTableCtrl<-gainCountsPerPatientTable[which(gainCountsPerPatientTable$Sample %in% patCtrl),]

write.table(gainCountsPerPatientTableCase, file="gainCountsPerPatientTableCase.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(gainCountsPerPatientTableCtrl, file="gainCountsPerPatientTableCtrl.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesCase<-as.matrix(apply(gainCountsPerPatientTableCase[, 2:dim(gainCountsPerPatientTableCase)[2]], 2, sum))
genesCtrl<-as.matrix(apply(gainCountsPerPatientTableCtrl[, 2:dim(gainCountsPerPatientTableCtrl)[2]], 2, sum))

write.table(genesCase, file="gainGenesCase.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesCtrl, file="gainGenesCtrl.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)

### loss

patCase<-rep("", dim(selCase)[1])
patCtrl<-rep("", dim(selCtrl)[1])

patList<-makelist(lossCountsPerPatientTable)

patCase<-unlist(patList[1])
patCtrl<-unlist(patList[2])

patCase<-patCase[which(patCase!="")]
patCtrl<-patCtrl[which(patCtrl!="")]

lossCountsPerPatientTableCase<-lossCountsPerPatientTable[which(lossCountsPerPatientTable$Sample %in% patCase),]
lossCountsPerPatientTableCtrl<-lossCountsPerPatientTable[which(lossCountsPerPatientTable$Sample %in% patCtrl),]

write.table(lossCountsPerPatientTableCase, file="lossCountsPerPatientTableCase.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
write.table(lossCountsPerPatientTableCtrl, file="lossCountsPerPatientTableCtrl.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

genesCase<-as.matrix(apply(lossCountsPerPatientTableCase[, 2:dim(lossCountsPerPatientTableCase)[2]], 2, sum))
genesCtrl<-as.matrix(apply(lossCountsPerPatientTableCtrl[, 2:dim(lossCountsPerPatientTableCtrl)[2]], 2, sum))

write.table(genesCase, file="lossGenesCase.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)
write.table(genesCtrl, file="lossGenesCtrl.tsv", sep="\t", dec=".", col.names=NA, row.names=T, quote=F)


