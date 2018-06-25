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

# DESCRIPTION

# THIS VERSION EXCLUDES ALL CUSTOM FILTERS. 
# VARSCAN OUTPUT .snp.pairs.tsv WERE FILTERD TO RETAIN ONLY SOMATIC MUTATIONS
# somaticFilter PROVIDED IN THE VARSCAN JAVA PACKAGE WILL BE APPLIED TO .tsv.somatic 

# somaticFilter
# This command filters somatic mutation calls to remove clusters of false positives and SNV calls near indels. Note: this is a basic filter. More advanced filtering strategies consider mapping quality, read mismatches, soft-trimming, and other factors when deciding whether or not to filter a variant. See the VarScan 2 publication (Koboldt et al, Genome Research, Feb 2012) for details.
# 
# USAGE: java -jar VarScan.jar somaticFilter [mutations file] OPTIONS
# mutations file - A file of SNVs from VarScan somatic
# 
# OPTIONS:
# --min-coverage  Minimum read depth [10]
# --min-reads2    Minimum supporting reads for a variant [2]
# --min-strands2  Minimum # of strands on which variant observed (1 or 2) [1]
# --min-avg-qual  Minimum average base quality for variant-supporting reads [20]
# --min-var-freq  Minimum variant allele frequency threshold [0.20]
# --p-value       Default p-value threshold for calling variants [1e-01]
# --indel-file    File of indels for filtering nearby SNPs
# --output-file   Optional output file for filtered variants

#library(BSgenome.Hsapiens.UCSC.hg19)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

# getRefBases<-function(annoSet){
#   
#   Chroms<-annoSet$Chr
#   sel<-which(Chroms=="23")
#   if(length(sel)>0) Chroms[sel]<-"X"
#   
#   sel<-which(Chroms=="24")
#   if(length(sel)>0) Chroms[sel]<-"Y"
#   
#   Chroms<-paste("chr", Chroms, sep="")
#   Starts<-annoSet$Start
#   Starts<-Starts
#   Ends<-Starts
#     
#   refs<-getSeq(Hsapiens, Chroms, start=Starts, end=Ends, as.character=T)
#   
#   return(refs)
# }

setEnds<-function(dset){
  
  #dset$Start<-dset$Start+1
  #Ends<-dset$Start
  Ends<-dset$End
  altSeq<-dset$Alt
  sel<-which(substr(altSeq, 1, 1)=="-")
  lenghts<-nchar(altSeq[sel])-2
  altSeq<-substr(altSeq, 2, nchar(altSeq))
  dset[sel, 'Alt']<-"-"
  dset[sel, 'Ref']<-altSeq[sel]
  dset[-sel, 'Alt']<-altSeq[-sel]
  dset[-sel, 'Ref']<-"-"
  Ends[sel]<-Ends[sel]+lenghts
  dset$End<-Ends
  
  return(dset)
}


convertVarscan2Annovar <- function () {
  tmpSet=NULL; annoSet=NULL;
  
  annovarDIR="~/annovar/"

  VARSCAN="~/tools/"

  #----------------------------------------------------------------------
  #----------------------------------------------------------------------
  
  # PROCESSING ARGUMENTS
  
  args <- commandArgs(trailingOnly = TRUE)
  
  arg<-args[1]
  
  #arg<-"/media/marco/DATA2/NGS/WES/Projects/2016_01_29/Results"
  
  setwd(arg)
  
  samplesFileName="samtoolsList"

  exGenesList<-"~/reference/genes_exclusion_list.csv"

  fileNames<-list.files(arg, pattern="\\.indels.pairs.tsv")
  
  filesNum<-length(fileNames)
  
  # LOADING AND PROCESSING DATA
  
  samplesList<-read.delim(paste(getwd(), samplesFileName, sep="/"), header=F, quote="", stringsAsFactors=F)
#############################################eugenio: changed header to FALSE###################################################

  exGenesList<-read.delim(exGenesList, header=F, quote="", stringsAsFactors=F)
  
  for(i in 1:filesNum){
    
    dset<-read.delim(paste(getwd(), fileNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    
    sel<-which((dset$somatic_status=="Somatic") | (dset$somatic_status=="LOH"))
    dset<-dset[sel,]
    
    write.table(dset, file=paste0(fileNames[i], ".somatic"), sep="\t", row.names=F, col.names=T, quote=F)
    
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  fileNames<-list.files(arg, pattern="\\.indels.pairs.tsv.somatic$")
  
  filesNum<-length(fileNames)
  
  for(i in 1:filesNum){
    
    dset<-read.delim(paste(getwd(), fileNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    
#     ### CALL VARSCAN somaticFilter
    
    CMDString<-paste0("java -jar ",  VARSCAN, "VarScan.v2.3.9.jar somaticFilter ", paste(getwd(), fileNames[i], sep="/"), 
                      " --min-coverage 1 ", 
                      " --min-reads2 2 ", 
                      " --min-var-freq 0.1 ", 
                      " --output-file ", paste0(getwd(), "/", fileNames[i], ".filtered"))
    
    system(CMDString)
    
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  fileNames<-list.files(arg, pattern="\\.indels.pairs.tsv.somatic.filtered$")
  
  filesNum<-length(fileNames)
  
    for(i in 1:filesNum){
    
    dset<-read.delim(paste(getwd(), fileNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    
    sel<-which(dset$chrom=="X")
    if(length(sel)>0) dset[sel, 'chrom']=23
    
    sel<-which(dset$chrom=="Y")
    if(length(sel)>0) dset[sel, 'chrom']=24
    
    ### Check two strands supporting reads
    
    sel<-which(dset$tumor_reads2_plus==0 | dset$tumor_reads2_minus==0)
    if(length(sel>0)) dset<-dset[-sel,]
    
    ### MAP TO SAMPLE NAME
    
    m<-dim(dset)[1]
    
    sampleName<-fileNames[i]
    #keyNum<-paste0(unlist(strsplit(unlist(strsplit(sampleName, "[.]")[1]), "_"))[3], "_", unlist(strsplit(unlist(strsplit(sampleName, "[.]")[1]), "_"))[4])
    #keyNum<-unlist(strsplit(sampleName, "_"))[4]
    #sel<-which(samplesList[1]==keyNum)
    
    #######eugenio: added following conditions, to distinguish between sample names like "Sample_1258_sequence_1.fastq" and sample names like "Sample_3FK_3D_sequence_1.fastq"
    keyNum<-unlist(strsplit(sampleName, "_"))
    print(length(keyNum)) ###############################################
    if(length(keyNum) == 11){
      keyNum = keyNum[4]
    }
    if(length(keyNum) == 13){
      keyNum = keyNum[5]
    }
    print(keyNum) ###############################################
    sel<-which(grepl(keyNum, samplesList[,1]))
    sampleLabel<-samplesList[sel,]
    
    labels<-rep(sampleLabel, m)
    
    ### BUILD CURRENT DATASET
    
    tmpSet<-cbind(dset$chrom, dset$position, dset$position, dset[, 3:23], labels, stringsAsFactors=F)
    namesVector<-colnames(tmpSet)
    namesVector[1]<-"Chr"
    namesVector[2]<-"Start"
    namesVector[3]<-"End"
    namesVector[4]<-"Ref"
    namesVector[5]<-"Alt"
    namesVector[25]<-"Sample"
    colnames(tmpSet)<-namesVector
    
    #annoSet$Ref<-getRefBases(annoSet)
    
    annoSet<-rbind(annoSet, tmpSet)
    
  }
########################eugenio: added an if condition, to prevent code from halting in case no calls met requirements
  if(dim(annoSet)[1] != 0){
    annoSet<-setEnds(annoSet)
  }
  
  ### WRITE OUT REPORT
  
  write.table(annoSet, file="annoSetVarscan_indels.tsv", row.names=F, col.names=F, sep="\t", quote=F)
  
  ### CALL ANNOVAR
  
  CMDString<-paste("perl ", annovarDIR, "table_annovar.pl ", arg, "/annoSetVarscan_indels.tsv ",  annovarDIR, "humandb/ -buildver hg19 -out ", arg, "/somatic_indels_varscan ", 
                   "-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,snp138,1000g2014oct_all,exac03nontcga,ljb26_all -operation g,r,r,f,f,f,f,f -nastring NA -csvout", sep="")
  
  system(CMDString)
  
   #CMDString<-paste("rm ", arg, "/annoSet.tsv", sep="")
   #system(CMDString)
  
  ### BUILD ANNOTATED DATASET
  
  annotatedSetFileName<-"somatic_indels_varscan.hg19_multianno.csv"
  
  tmpSet<-read.delim(paste(getwd(), annotatedSetFileName, sep="/"), header=T, sep=",", stringsAsFactors=F)
  
  outSet<-cbind(annoSet$Sample, tmpSet[,1:16], tmpSet[,c(24:25, 27, 35)], annoSet[,6:14], annoSet[,15:24])
  namesVector<-colnames(outSet)
  namesVector[1]<-"Sample"
  colnames(outSet)<-namesVector
  
  ### REMOVE NON-EXOMIC, SYNONIMOUS
  
  sel<-which((outSet$Func.refGene=="exonic") | (outSet$Func.refGene=="exonic;exonic"))
  outSet<-outSet[sel,]
  sel<-which((outSet$ExonicFunc.refGene=="unknown;unknown") | (outSet$ExonicFunc.refGene=="unknown"))
  #             (outSet$ExonicFunc.refGene=="unknown") | (outSet$ExonicFunc.refGene==""))
  if(length(sel>0)) {outSet<-outSet[-sel,]}
  
  ### REMOVE EXCLUSION LIST GENES
  
  pattern<-paste("^", exGenesList[,1], sep="")
  for(i in 1:length(pattern)){
    sel<-grep(pattern[i], outSet$Gene.refGene)
    if(length(sel)>0) {
      write.table(outSet[sel,], file="exclusion_list_indels_removed.tsv", row.names=F, col.names=T, sep="\t", dec=".", quote=F)
      outSet<-outSet[-sel,]
    }
  }
  
  
  write.table(outSet, file="somatic_indels_varscan_hg19_annotated.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")
  
  CMDString<-paste("rm ", arg, "/", annotatedSetFileName, sep="")
  
  system(CMDString)
  
  CMDString<-paste0("rm ", arg, "/*.somatic")
  
  system(CMDString)
  
  CMDString<-paste0("rm ", arg, "/*.somatic.filtered")
  
  system(CMDString)
}

convertVarscan2Annovar()
