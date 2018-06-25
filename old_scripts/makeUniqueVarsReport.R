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

# PROCESSING ARGUMENTS

args <- commandArgs(trailingOnly = TRUE)

arg<-args[1]

#arg<-"/media/marco/DATA2/NGS/WES/Projects/2016_01_29/Results"

setwd(arg)

varscanSNVFileName<-"somatic_mutation_varscan_hg19_annotated.tsv"
mutectSNVFileName<-"somatic_mutation_hg19_annotated.tsv"
varscanINDELSFileName<-"somatic_indels_varscan_hg19_annotated.tsv"

makeUniqueVarsReport <- function () {
  varscanSNVFile<-read.delim(paste(getwd(), varscanSNVFileName, sep="/"), header=T, quote="", stringsAsFactors=F)
  mutectSNVFile<-read.delim(paste(getwd(), mutectSNVFileName, sep="/"), header=T, quote="", stringsAsFactors=F)
  varscanINDELSFile<-read.delim(paste(getwd(), varscanINDELSFileName, sep="/"), header=T, quote="", stringsAsFactors=F)
  
  # ADD COLUMN WITH DETECTION ALGORITHM
  
  varscanMethod<-rep("VARSCAN", dim(varscanSNVFile)[1])
  varscanSNVFile<-cbind(varscanSNVFile, varscanMethod, stringsAsFactors=F)
  colnames(varscanSNVFile)[length(colnames(varscanSNVFile))]<-"DETECTION_TYPE"
  
  varscanMethod<-rep("VARSCAN", dim(varscanINDELSFile)[1])
  varscanINDELSFile<-cbind(varscanINDELSFile, varscanMethod, stringsAsFactors=F)
  colnames(varscanINDELSFile)[length(colnames(varscanINDELSFile))]<-"DETECTION_TYPE"
  
  mutectMethod<-rep("MUTECT", dim(mutectSNVFile)[1])
  mutectSNVFile<-cbind(mutectSNVFile, mutectMethod, stringsAsFactors=F)
  colnames(mutectSNVFile)[length(colnames(mutectSNVFile))]<-"DETECTION_TYPE"
  
  # PREPARING FILES FOR MERGE
  
  namesVector<-c("SAMPLE", "CHR", "START", "END", "REF", "ALT", "EXON", "GENE", "SOMATIC_STATUS", "VAR_TYPE", "AA_CHANGE", "CYTOBAND", "SEGDUPS","ESP", 
                 "DBSNP_ID", "TGP", "EXAC", "SIFT_SCORE", "SIFT_PRED", "POLY_PRED", "MUTASSES_PRED",
                 "N_FRACTION", "T_FRACTION", "N_REF", "N_ALT", "T_REF", "T_ALT", "DETECTION_METHOD")
  
  # VARSCAN SNV
  
  normal_fraction<-round(varscanSNVFile$normal_reads2/(varscanSNVFile$normal_reads1+varscanSNVFile$normal_reads2), digits=2)
  
  tumor_fraction<-round(varscanSNVFile$tumor_reads2/(varscanSNVFile$tumor_reads1+varscanSNVFile$tumor_reads2), digits=2)
  
  tmpVarscanSNVFile<-cbind(varscanSNVFile[, 1:8], varscanSNVFile[, 'somatic_status'], varscanSNVFile[, 10:21], normal_fraction, tumor_fraction, varscanSNVFile$normal_reads1, varscanSNVFile$normal_reads2, varscanSNVFile$tumor_reads1, varscanSNVFile$tumor_reads2, varscanSNVFile[, 'DETECTION_TYPE'], stringsAsFactors=F)
  
  colnames(tmpVarscanSNVFile)<-namesVector
  
  # MUTECT SNV
  
  normal_fraction<-round(mutectSNVFile$n_alt_count/(mutectSNVFile$n_ref_count+mutectSNVFile$n_alt_count), digits=2)
  
  tumor_fraction<-round(mutectSNVFile$t_alt_count/(mutectSNVFile$t_ref_count+mutectSNVFile$t_alt_count), digits=2)
  
  tmpMutectSNVFile<-cbind(mutectSNVFile[, 1:8], mutectSNVFile[, 'somatic_status'], mutectSNVFile[, 10:21], normal_fraction, tumor_fraction, mutectSNVFile$n_ref_count, mutectSNVFile$n_alt_count, mutectSNVFile$t_ref_count, mutectSNVFile$t_alt_count, mutectSNVFile[, 'DETECTION_TYPE'], stringsAsFactors=F)
  
  colnames(tmpMutectSNVFile)<-namesVector
  
  # VARSCAN INDELS
  
  normal_fraction<-round(varscanINDELSFile$normal_reads2/(varscanINDELSFile$normal_reads1+varscanINDELSFile$normal_reads2), digits=2)
  
  tumor_fraction<-round(varscanINDELSFile$tumor_reads2/(varscanINDELSFile$tumor_reads1+varscanINDELSFile$tumor_reads2), digits=2)
  
  tmpVarscanINDELSFile<-cbind(varscanINDELSFile[, 1:8], varscanINDELSFile[, 'somatic_status'], varscanINDELSFile[, 10:21], normal_fraction, tumor_fraction, varscanINDELSFile$normal_reads1, varscanINDELSFile$normal_reads2, varscanINDELSFile$tumor_reads1, varscanINDELSFile$tumor_reads2, varscanINDELSFile[, 'DETECTION_TYPE'], stringsAsFactors=F)
  
  colnames(tmpVarscanINDELSFile)<-namesVector
  
  # MERGING FILE AND SORTING BY SAMPLE
  
  tmpUniqueVariants=NULL
  
  tmpUniqueVariants<-rbind(tmpMutectSNVFile, tmpVarscanSNVFile, tmpVarscanINDELSFile)
  
  tmpUniqueVariants<-tmpUniqueVariants[order(tmpUniqueVariants$SAMPLE, decreasing=T),]
  
  # FIND DUPLICATES
  
  dup<-duplicated(tmpUniqueVariants[, 1:6])
  
  dupValidated<-which(dup==TRUE)
  
  if(length(dupValidated)>0){
  
    for(i in 1:length(dupValidated)){
      
      sel<-which(tmpUniqueVariants$GENE==tmpUniqueVariants[dupValidated[i], 'GENE'] & tmpUniqueVariants$SAMPLE==tmpUniqueVariants[dupValidated[i], 'SAMPLE'] &  
      tmpUniqueVariants$START==tmpUniqueVariants[dupValidated[i], 'START'] & tmpUniqueVariants$END==tmpUniqueVariants[dupValidated[i], 'END'] & tmpUniqueVariants$REF==tmpUniqueVariants[dupValidated[i], 'REF'] &
      tmpUniqueVariants$ALT==tmpUniqueVariants[dupValidated[i], 'ALT'])
    
      tmpUniqueVariants[sel[1], 'DETECTION_METHOD']="BOTH"
    
    }
  
  }
  tmpUniqueVariants[is.na(tmpUniqueVariants)]=""
  uniqueVariants<-tmpUniqueVariants[!dup,]
  
  # FIND PER-GENE REPEATED MUTATIONS
  rep=data.frame()
  distances=NULL
  for(j in 1:length(levels(as.factor(uniqueVariants$SAMPLE)))){
    
    duniqueVariants<-uniqueVariants[uniqueVariants$SAMPLE==levels(as.factor(uniqueVariants$SAMPLE))[j],]
    
    dupGeneVariants<-duplicated(duniqueVariants[, c(1, 8)])
    selDup<-duniqueVariants[dupGeneVariants,]
    for(i in 1:length(levels(as.factor(selDup$GENE)))){
    
      sel<-which(duniqueVariants$GENE==levels(as.factor(selDup$GENE))[i])
      if(length(sel)>0){
        if(length(sel)[1]==2){
      
          d=abs(duniqueVariants[sel[2], 'START'] - duniqueVariants[sel[1], 'START'])
        
          cat(duniqueVariants[sel[1], 'GENE'], "\n")
          cat(d, "\n")
          
          distances<-c(distances, rep(d, length(sel)))
          rep=rbind(rep, duniqueVariants[sel,])
      
        } else {
          
          d=0
          distances<-c(distances, d)
          
          l=length(sel)
          for(g in 2:l) {
          #d=abs(ceiling((duniqueVariants[sel[length(sel)], 'START'] - duniqueVariants[sel[1], 'START'])/length(sel)))
            d=abs(duniqueVariants[sel[g], 'START'] - duniqueVariants[sel[g-1], 'START'])
          
            cat(duniqueVariants[sel[1], 'GENE'], "\n")
            cat(d, "\n")
          
            distances<-c(distances, d)
            }
          rep=rbind(rep, duniqueVariants[sel,])
        }
      }
    }
  }
  
  distances<-data.frame(distances, rep$GENE)

  ### removing hypermutated

  sel<-which(distances[,1]>0 & distances[,1]<=100)
  if(length(sel)>0){
    
    for(i in 1:length(sel)){
      
      remove<-which(uniqueVariants$SAMPLE==rep[sel[i], 'SAMPLE'] & uniqueVariants$CHR==rep[sel[i], 'CHR'] & 
        uniqueVariants$START==rep[sel[i], 'START'] & uniqueVariants$END==rep[sel[i], 'END'])
        # & uniqueVariants$REF==rep[sel[i], 'REF'] & uniqueVariants$ALT==rep[sel[i], 'ALT'])
    
      #uniqueVariants<-uniqueVariants[-remove,]
      
    }
    
    ##
    
    write.table(rep[sel,], file="hyper_removed.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=T)
    
    rep<-rep[-sel,]
  }
  
  
    
  # removing high frequency mutated genes rsID associated
  # <=2 mutations per-gene both rsID associated
  # >2 mutations per-gene and 50% rsID associated  
  
  #polymorph<-rep(0, dim(rep)[1])
  
  #polymorph<-rep(0, dim(uniqueVariants)[1])
  
  #polyrem<-data.frame()
    
  #for(i in 1:length(levels(as.factor(rep$GENE)))){
  #for(i in 1:length(levels(as.factor(uniqueVariants$GENE)))){
      
      #sel<-which(rep$GENE==levels(as.factor(rep$GENE))[i])
  #    sel<-which(uniqueVariants$GENE==levels(as.factor(uniqueVariants$GENE))[i])
  #    if(length(sel)>1){
  #      if(length(sel)[1]==2){
          
          #nsnp=length(which(rep[sel, 'DBSNP_ID']!="."))
  #        nsnp=length(which(uniqueVariants[sel, 'DBSNP_ID']!="."))
  #        if(nsnp==2) {polymorph[sel]=TRUE}
                    
  #      } else {
          
          #nsnp=length(which(rep[sel, 'DBSNP_ID']!="."))
  #        nsnp=length(which(uniqueVariants[sel, 'DBSNP_ID']!="."))
  #        if(nsnp/length(sel)>=0.5) {polymorph[sel]=TRUE}
          
  #      }
  #    }
  #  }
  
  #sel<-which(polymorph==TRUE)
  #if(length(sel)>0){
  
#     for(i in 1:length(sel)){
#         remove<-which(uniqueVariants$SAMPLE==rep[sel[i], 'SAMPLE'] & uniqueVariants$CHR==rep[sel[i], 'CHR'] & 
#                     uniqueVariants$START==rep[sel[i], 'START'] & uniqueVariants$END==rep[sel[i], 'END'] &
#                     uniqueVariants$REF==rep[sel[i], 'REF'] & uniqueVariants$ALT==rep[sel[i], 'ALT'])
#         
#         if(length(remove)>0) {  polyrem<-rbind(polyrem, uniqueVariants[remove,])  
#                                 uniqueVariants<-uniqueVariants[-remove,]
#         }
#       
#       
#       }
    
   # polyrem<-uniqueVariants[sel,]
    
   # uniqueVariants<-uniqueVariants[-sel,]
    
   # }
  
  #if(dim(polyrem)[1]>0) {write.table(polyrem, file="polymorphic_removed.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=T)}
  if(dim(uniqueVariants)[1]>0) {write.table(uniqueVariants, file="somatic_variants.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")}
  if(dim(rep)[1]>0) {write.table(rep, file="somatic_variants_repeated.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")}
  if(dim(distances)[1]>0) {write.table(distances, file="dist.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")}
  
}

makeUniqueVarsReport()
