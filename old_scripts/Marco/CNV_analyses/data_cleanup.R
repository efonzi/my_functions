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

# FUNCTIONS

makeDgvFile<-function(dgvFile, rev){
  
  dgv<-read.delim(paste(DGVFILENAME, ".tsv", sep=""), header=T, sep=" ", stringsAsFactors=F)
  dgv<-dgv[which(dgv$varianttype=="CNV"),]
  dgv[which(dgv$chr=="X"), "chr"]<-23
  dgv[which(dgv$chr=="Y"), "chr"]<-24
  dgv<-dgv[which(dgv$chr!="5_h2_hap1"),]
  dgv<-dgv[which(dgv$chr!="6_cox_hap1"),]  
  
  write.table(file=paste(DGVFILENAME, "_", rev, ".tsv", sep=""), dgv, sep="\t", row.names=F, col.names=T)
  
}

createRepTable<-function(repTab){
  
  chr<-matrix(unlist(strsplit(as.character(repTab[,1]), split=c(":"), fixed=T)), ncol=2, byrow=T)
  
  chr[,1]<-gsub("chr", "", chr[,1])
  
  chr[which(chr[, 1]=="X"), 1]<-23
  
  chr[which(chr[, 1]=="Y"), 1]<-24
  
  pos<-matrix(unlist(strsplit(as.character(chr[,2]), split=c("-"), fixed=T)), ncol=2, byrow=T)
  
  pos[,1]<-gsub(",", "", pos[, 1])
  
  pos[,2]<-gsub(",", "", pos[, 2])
  
  cnv<-repTab[,2]
  
  cnv[cnv=="CN Gain"]<-"gain"
  
  cnv[cnv=="CN Loss"]<-"loss"
  
  cnv[cnv=="High Copy Gain"]<-"duplication"
  
  cnv[cnv=="Homozygous Copy Loss"]<-"deletion"
  
  repTab<-data.frame(as.numeric(chr[,1]), as.numeric(pos[,1]), as.numeric(pos[,2]), repTab$Probe.Median, cnv, repTab$Sample, repTab$Gene.Symbol, stringsAsFactors=F)
  
  colnames(repTab)<-c("CHR", "START", "END", "MEDIAN", "CNV_TYPE", "SAMPLE", "GENE")
  
  return(repTab)
}

createDGVDeletedSetDebug<-function(){
  
  #FOR DEBUG PURPOSE A FULL DATASET WILL BE CREATED
  #IT REQUIRE MORE DISK SPACE (200MB WITH AN INPUT REPORT OF 36000 ENTRIES
  
  if(dim(sel)[1]!=0){
    #found=0
    for(j in 1:dim(sel)[1]){
    
      if (sel[j, "variantsubtype"]==tab[i, "CNV_TYPE"]) {
      
        remove[i]="TRUE"
        repDeleted<-data.frame(NULL)
        for(j in 1:dim(sel)[1]){
          repDeleted<-rbind(repDeleted, tab[i,])
        }
        repDeleted<-cbind(repDeleted, sel)
        dgvDeleted<-rbind(dgvDeleted, repDeleted)
        cat(dim(dgvDeleted)[1], "\n")
        break
      }
    }
  }
  return(list(dgvDeleted=dgvDeleted, remove=remove))
}

createDGVDeletedSet<-function(sel, i, remove, dgvDeleted){

  if(dim(sel)[1]!=0){
    #found=0
    for(j in 1:dim(sel)[1]){
    
      if (sel[j, "variantsubtype"]==tab[i, "CNV_TYPE"]) {
      
        remove[i]="TRUE"
        dgvDeleted<-rbind(dgvDeleted, tab[i, ])
        #cat(dim(dgvDeleted))
        break
      }
    }
  }
  return(list(dgvDeleted=dgvDeleted, remove=remove))
}

mergeNexusReports<-function(REP1, REP2){
  
  TAB1<-read.delim(paste(REP1, ".txt", sep="") , header=T, sep="\t", stringsAsFactors=F)
  TAB2<-read.delim(paste(REP2, ".txt", sep="") , header=T, sep="\t", stringsAsFactors=F)
  REP<-rbind(TAB1, TAB2)
  return(REP)
}

filter_set<-function(dset){
  sel<-which(dset$Event=="Allelic Imbalance")
  if(length(sel)>0){
    dset<-dset[-sel, ]
  }
  sel<-which(dset$Probes<5)
  if(length(sel)>0){
    dset<-dset[-sel,]
  }
  sel<-which(dset$Length<1000)
  if(length(sel)>0){
    dset<-dset[-sel,]
  }
#   sel<-which((dset$Event=="gain") & (dset$Length<1000))
#   if(length(sel)>0){
#     dset<-dset[-sel, ]
#   }
#   sel<-which((dset$Event=="duplication") & (dset$Length<1000))
#   if(length(sel)>0){
#     dset<-dset[-sel, ]
#   }
#   sel<-which((dset$Event=="loss") & (dset$Length<1000))
#   if(length(sel)>0){
#     dset<-dset[-sel, ]
#   }
#   sel<-which((dset$Event=="deletion") & (dset$Length<1000))
#   if(length(sel)>0){
#     dset<-dset[-sel, ]
#   }
#   sel<-which((dset$Event=="LOH") & (dset$Length<1000))
#   if(length(sel)>0){
#     dset<-dset[-sel, ]
#   }
  return(dset)
}

filter <- function (tab) {
  
  dgvDeleted<-data.frame(NULL)
  
  remove<-rep("FALSE", dim(tab)[1])
  
  for(i in 1:dim(tab)[1]){
    
    cat("Checking event #", i, " ", tab[i, "Chromosome.Region"], "\n")
    
    sel<-dgv[which(dgv$chr==tab[i, "CHR"]),]
    
    # SELECTION BASED ON START & END POSITIONS (ALL REGIONS INSIDE THE INTERVAL WILL BE SELECTED AND FILTERED OUT)
    #-------------------------------------------------------------------------------------------------------------
    
    sel<-sel[which((sel$start<=tab[i, "START"]-s1) & (sel$end>=tab[i, "END"]+e2) & (sel$variantsubtype==tab[i, 'CNV_TYPE'])),] 
    
    
    # SELECTIONS BASED ON START & END POSITION (# BP TOLERANCE) ONLY LENGTH-MATCHING REGIONS WILL BE SELECTED AND FILTERED OUT
    #---------------------------------------------------------------------------------------------------------------------------
    #   
    #  sel<-sel[which((sel$start>=tab[i, "START"]-s1) & (sel$start<=tab[i, "START"]+s2)),]
    #   
    #  sel<-sel[which((sel$end>=tab[i, "END"]-e1) & (sel$end<=tab[i, "END"]+e2)),]
    
    
    
    delList<-createDGVDeletedSet(sel, i, remove, dgvDeleted)
    
    dgvDeleted<-delList$dgvDeleted
    remove<-delList$remove
    
  }
  
  cat(paste("Entries to be removed: ", length(which(remove=="TRUE")), sep=""), "\n")
  
  write.table(dgvDeleted, file=paste(REPFILENAME, "_dgv_deleted.tsv", sep=""), sep="\t", col.names=T, row.names=T)
  
  # DELETE CNV FROM REPORT FILE
  
  repTab<-repTab[remove=="FALSE",]
  
  write.table(repTab, file=paste(REPFILENAME, "_dgv_filtered.tsv", sep=""), sep="\t", col.names=T, row.names=F)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#setwd("/home/mmanfrini/data1/data/Bologna/ProjectAE")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/chromothripsis")

#REPFILENAME="snp-recovery-list.tsv"
REPFILENAME="missing_sample.tsv"

#DGVPATH="/media/DATA1/DATA/Genomes/DGV"
DGVPATH="/Users/Eugenio/Desktop/Lavoro/Seragnoli/annotation_files/DGV"
DGVFILENAME="var_hg19_1.tsv"
NREFFILENAME="normal_ref_HD.tsv"

# LOAD REPORT FILE

repTab<-read.delim(paste(REPFILENAME, sep="/") , header=T, sep="\t", stringsAsFactors=F)

#repTab<-mergeNexusReports("report APL 6.0", "REPORT APL CYTOSCAN")

# LOAD DGV  & REF FILE

dgv<-read.delim(paste(DGVPATH, DGVFILENAME, sep="/") , header=T, sep="", stringsAsFactors=F)
dgv<-dgv[, c(2, 3, 4 ,6)]
nrefhd<-read.delim(paste(DGVPATH, NREFFILENAME, sep="/") , header=T, sep="", stringsAsFactors=F)

dgv<-rbind(dgv, nrefhd)

# FILTER DATASET

repTab<-filter_set(repTab)

# CREATE TABLE FROM REPORT FILE WITH CHR, START, END, CNV TYPE INFOS

tab<-createRepTable(repTab)

#tab<-tab[1:1000,] #FOR TESTING PURPOSE ONLY

# CREATE DATA FRAMES TO STORE INFOS FROM DELETED ENTRIES

# SEARCH FOR CNV IN DGV REF FILE (EXACT LENGTH-MATCHING)

#      s   S   s                             e   E   e
# DGV  |---|---|-----------------------------|---|---|
# REP      |-------------------------------------|
#
#
# S±s=±100bp, E±e=±100bp;


# SEARCH FOR CNV IN DGV REF FILE (FLANKING REGIONS EXTENSION)

#      s1   S   s2                            e1  E   e2
# DGV  |----|---|-----------------------------|---|----|
# REP       |-------------------------------------|
#
#
# s1-S=-1000bp, s2-S=1000bp; e1-E=-1000bp; e2-E=1000bp

# # BP TOLERANCE

s1=1000
s2=1000
e1=1000
e2=1000

# USUAL INTERVALS - TESTED OK
# 
# s1=50000
# s2=50000
# e1=50000
# e2=50000

#-----------------------------------
# SUBTRACT DGV & REF
#-----------------------------------

filter(tab)

