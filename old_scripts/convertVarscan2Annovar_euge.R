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


#-------------------------------------------------------------------------------------------------------
#---------------------------------------------- FUNCTIONS ----------------------------------------------
#-------------------------------------------------------------------------------------------------------


# adjust "ref" and "alt" for indels, in order to be compatible with ANNOVAR desired input
setEnds<-function(dset){
    
    Ends<-dset$End
    altSeq<-dset$Alt
    
    #select ALT that start by "-", which are the DELETIONS
    sel<-which(substr(altSeq, 1, 1)=="-")
    #save the lengths of all DELETIONS minus 2
    lenghts<-nchar(altSeq[sel])-2
    #remove first character from all ALT (which is either "+" or "-")
    altSeq<-substr(altSeq, 2, nchar(altSeq))
    
#    cat("1\n")
    
    #in ALT column of dset, substitute all deletions with "-"
    dset[sel, 'Alt']<-"-"
    
    #in REF column of dset, substitute all deletions with deleted string (without "-")
    dset[sel, 'Ref']<-altSeq[sel]
    
#    cat("2\n")
    
    if(length(sel)>0){
      
      #in ALT column of dset, substitute all insertions with inserted string (without "+")
      dset[-sel, 'Alt']<-altSeq[-sel]
      #in REF column of dset, substitute all insertions with "-"
      dset[-sel, 'Ref']<-"-"
    }
    else{ # in case of no deletions, "sel" will be empty and the following code is needed
      
      #in ALT column of dset, substitute all insertions with inserted string (without "+")
      dset[, 'Alt'] = altSeq
      #in REF column of dset, substitute all insertions with "-"
      dset[, 'Ref'] = "-"
    }
    
    #adjust END positions of deletions
    Ends[sel]<-Ends[sel]+lenghts
    dset$End<-Ends
    
    return(dset)
}


convertVarscan2Annovar <- function () {
  
  # directories and files
  annovarDIR="/media/DATA1/DATA/Software/annovar/"
  VARSCAN="/media/DATA1/DATA/Software/"
  exGenesList<-"/media/DATA1/DATA/Genomes/genes_exclusion_list.csv"
  #VARSCAN="/Users/Eugenio/Desktop/NGScli/varscan2.3.9/"
  #setwd("/Volumes/UNTITLED/varscan_out_temp")
  
  # PROCESSING ARGUMENTS
  
  args <- commandArgs(trailingOnly = TRUE)
  arg<-args[1]
  setwd(arg)



  #--------------------------------------------------------------------------------------------------------
  
  # LOAD DATA, select "somatic" and "loh" variants
  indelNames<-list.files(arg, ".indels.pairs$")
  snpNames<-list.files(arg, ".snp.pairs$")
  
  for(i in 1:length(snpNames)){
    
    snpSet<-read.delim(paste(getwd(), snpNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    indelSet<-read.delim(paste(getwd(), indelNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    
    snpSet<- snpSet[which((snpSet$somatic_status=="Somatic") | (snpSet$somatic_status=="LOH")),]
    indelSet<- indelSet[which((indelSet$somatic_status=="Somatic") | (indelSet$somatic_status=="LOH")),]
    
    write.table(snpSet, file=paste0(snpNames[i], "_somatic.tsv"), row.names=F, col.names=T, sep="\t", quote=F)
    write.table(indelSet, file=paste0(indelNames[i], "_somatic.tsv"), row.names=F, col.names=T, sep="\t", quote=F)
    
  }
  
  
  # load "somatic" variants and run "somaticFilter" from Varscan2.3.9
  indelNames<-list.files(arg, ".indels.pairs_somatic.tsv$")
  snpNames<-list.files(arg, ".snp.pairs_somatic.tsv$")
  
  for(i in 1:length(snpNames)){
    
    snpString<-paste0("java -jar ",  VARSCAN, "VarScan.v2.3.9.jar somaticFilter ", paste(getwd(), snpNames[i], sep="/"),
                      " --indel-file ", paste(getwd(), indelNames[i], sep="/"), 
                      " --min-coverage 1 ", 
                      " --min-reads2 2 ", 
                      " --min-var-freq 0.1 ", 
                      " --output-file ", paste0(getwd(), "/", snpNames[i], ".filtered.tsv"))
    
    
    indelString<-paste0("java -jar ",  VARSCAN, "VarScan.v2.3.9.jar somaticFilter ", paste(getwd(), indelNames[i], sep="/"),
                        " --min-coverage 1 ",
                        " --min-reads2 2 ",
                        " --min-var-freq 0.1 ",
                        " --output-file ", paste0(getwd(), "/", indelNames[i], ".filtered.tsv"))
    
    
    system(snpString)
    system(indelString)
    
  }

  
  # load filtered files, change "X" and "Y" with "23" and "24",
  # remove all variants with "tumor_reads2" = 0, 
  # adjust "REF" and "ALT" to ANNOVAR standard (only for indels),
  # concatenate all results in one single data.frame called "annoSet"
  indelNames<-list.files(arg, "indels.pairs_somatic.tsv.filtered.tsv$")
  snpNames<-list.files(arg, "snp.pairs_somatic.tsv.filtered.tsv$")

  annoSet=NULL
  for(i in 1:length(snpNames)){
    
    snpSet<-read.delim(paste(getwd(), snpNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    cat("snp", i, "\n")
    indelSet<-read.delim(paste(getwd(), indelNames[i], sep="/"), header=T, quote="", stringsAsFactors=F)
    cat("indel", i, "\n")
    
    sel<-which(snpSet$chrom=="X")
    if(length(sel)>0) snpSet[sel, 'chrom']=23
    sel<-which(indelSet$chrom=="X")
    if(length(sel)>0) indelSet[sel, 'chrom']=23
    
    sel<-which(snpSet$chrom=="Y")
    if(length(sel)>0) snpSet[sel, 'chrom']=24
    sel<-which(indelSet$chrom=="Y")
    if(length(sel)>0) indelSet[sel, 'chrom']=24
    
    ### Check two strands supporting reads
    
    sel<-which(snpSet$tumor_reads2_plus==0 | snpSet$tumor_reads2_minus==0)
    if(length(sel>0)) snpSet<-snpSet[-sel,]
    sel<-which(indelSet$tumor_reads2_plus==0 | indelSet$tumor_reads2_minus==0)
    if(length(sel>0)) indelSet<-indelSet[-sel,]
    
    ### extract sample name from file name
    sampleName<-strsplit(snpNames[i], "[.]")[[1]][1]
    
    ### BUILD TEMPORARY DATASET FOR SNP

    labels<-rep(sampleName, dim(snpSet)[1])

    tmpSet<-cbind(snpSet$chrom, snpSet$position, snpSet$position, snpSet[, 3:23], labels, stringsAsFactors = FALSE)
    namesVector<-colnames(tmpSet)
    namesVector[1]<-"Chr"
    namesVector[2]<-"Start"
    namesVector[3]<-"End"
    namesVector[4]<-"Ref"
    namesVector[5]<-"Alt"
    namesVector[25]<-"Sample"
    colnames(tmpSet)<-namesVector

    annoSet<-rbind(annoSet, tmpSet)
    
    ### BUILD TEMPORARY DATASET FOR INDEL

    labels<-rep(sampleName, dim(indelSet)[1])

    tmpSet<-cbind(indelSet$chrom, indelSet$position, indelSet$position, indelSet[, 3:23], labels, stringsAsFactors = FALSE)
    namesVector<-colnames(tmpSet)
    namesVector[1]<-"Chr"
    namesVector[2]<-"Start"
    namesVector[3]<-"End"
    namesVector[4]<-"Ref"
    namesVector[5]<-"Alt"
    namesVector[25]<-"Sample"
    colnames(tmpSet)<-namesVector
    
    if(dim(tmpSet)[1]>0){
      cat("yes\n")
      tmpSet = setEnds(tmpSet) # apply function setEnds that adjust END, REF and ALT columns of INDEL variants
    }
    annoSet<-rbind(annoSet, tmpSet)
    cat(dim(annoSet)[1],"\n")
  }
  
  # some "T" are read as "TRUE" and this code fixes this
  for(r in 1:dim(annoSet)[1]){
    for(c in 1:dim(annoSet)[2]){
      if(annoSet[r,c] == "TRUE"){
        annoSet[r,c] = "T"}}}
  
  
  ### CREATE VARSCAN FOLDER AND WRITE OUT REPORT

  system(paste0("mkdir ", arg, "/varscan"))
  setwd(paste0(arg, "/varscan"))
  write.table(annoSet, file="annoSetVarscan.tsv", row.names=F, col.names=T, sep="\t", quote=F)

  
  
  #-------------------------------------- ANNOVAR CALL - DIFFERENT IN CASE OF MITOCHONDRIAL VARIANTS ----------------------------------------
  
  if(!"MT" %in% annoSet$Chr){
    
    ### CALL TABLE_ANNOVAR.PL
    
    CMDString<-paste("perl ", annovarDIR, "table_annovar.pl ", getwd(), "/annoSetVarscan.tsv ",  annovarDIR, "humandb/ -buildver hg19 -out ", arg, "/somatic_mutation_varscan ",
                     "-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,snp138,1000g2014oct_all,exac03nontcga,ljb26_all -operation g,r,r,f,f,f,f,f -nastring NA -csvout", sep="")
    
    system(CMDString)
    
    
    ### BUILD ANNOTATED DATASET
    
    annotatedSetFileName<-"somatic_mutation_varscan.hg19_multianno.csv"
    
    tmpSet<-read.delim(paste(getwd(), annotatedSetFileName, sep="/"), header=T, sep=",", stringsAsFactors=F)
    
    outSet<-cbind(annoSet[,25], tmpSet[,1:16], tmpSet[,c(24:25, 27, 35)], annoSet[,6:14], annoSet[,15:24])
    namesVector<-colnames(outSet)
    namesVector[1]<-"Sample"
    colnames(outSet)<-namesVector
    ### REMOVE NON-EXOIC, SYNONIMOUS
    
    sel<-which((outSet$Func.refGene=="exonic") | (outSet$Func.refGene=="exonic;exonic"))
    outSet<-outSet[sel,]
    sel<-which((outSet$ExonicFunc.refGene=="synonymous SNV") | (outSet$ExonicFunc.refGene=="synonymous SNV;synonymous SNV") |
                 (outSet$ExonicFunc.refGene=="unknown") | (outSet$ExonicFunc.refGene==""))
    outSet<-outSet[-sel,]
    
    
    ### REMOVE EXCLUSION LIST GENES
    
    exGenesList<-read.delim(exGenesList, header=F, quote="", stringsAsFactors=F)
    removed=data.frame()
    
    pattern<-paste("^", exGenesList[,1], sep="")
    for(i in 1:length(pattern)){
      sel<-grep(pattern[i], outSet$Gene.refGene)
      if(length(sel)>0) {
        removed<-rbind(removed, outSet[sel,])
        outSet<-outSet[-sel,]
      }
    }
    
    write.table(removed, file="exclusion_list_snv_removed.tsv", row.names=F, col.names=T, sep="\t", dec=".", quote=F)
    write.table(outSet, file="somatic_mutation_varscan_hg19_annotated.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")
    
    ############# REMOVE UNWANTED INTERMEDIATE FILES !!!!!!!!!! ###############
    
  }
  
  else{
    
    ### CALL ANNOTATE_VARIATION.PL FOR MITOCHONDRIAL VARIANTS
    
    CMDString<-paste0("perl ", annovarDIR, "annotate_variation.pl -buildver GRCh37_MT -dbtype ensGene ", getwd(), "/annoSetVarscan.tsv ", annovarDIR, "humandb/")
    system(CMDString)
    
    ### LOAD ANNOTATED EXONIC VARIANTS AND ADD COLUMN NAMES (ANNOVAR OUTPUT HAS NO NAMES)
    # write final table and delete intermediate files

    tmpSet<-read.delim("annoSetVarscan.tsv.exonic_variant_function", header=T, sep="\t", stringsAsFactors=F)
    
    outSet = tmpSet[,c(28,2:27)]
    colnames(outSet) = c("Sample", "variant_type", "annotation", colnames(annoSet)[1:24])
    
    write.table(outSet, file="mitochondrial_varscan_hg19_annotated.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")
    
    system(paste0("rm ", getwd(), "/annoSetVarscan.tsv.exonic_variant_function"))
    system(paste0("rm ", arg, "/*somatic*"))
    system(paste0("rm ", arg, "/*filtered*"))
    
  }
  
}

convertVarscan2Annovar()



