# Copyright (C) 2016  Marco Manfrini, PhD

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

library(rsnps)
Sys.setenv(http_proxy="http://antonella.padella2%40unibo.it:attila09.@192.168.65.4:8080")

# PROCESSING ARGUMENTS

args <- commandArgs(trailingOnly = TRUE)

arg<-args[1]

setwd(arg)

fileName="somatic_variants.tsv"

incGenesList<-"~/reference/genes_inclusion_list.csv"

incGenesList<-read.delim(incGenesList, header=F, quote="", stringsAsFactors=F)

# FILTER AGAINST ESP, TGP AND EXAC VAF

uniqueVariants<-read.delim(file=paste0(getwd(), "/", fileName), sep="\t", dec=".", stringsAsFactors=F)

print(paste("uno", dim(uniqueVariants)[1], sep = "_"))

sel<-which(uniqueVariants$EXAC>=0.01)
if(length(sel)>0){uniqueVariants<-uniqueVariants[-sel,]}

sel<-which(uniqueVariants$ESP>=0.01)
if(length(sel)>0){uniqueVariants<-uniqueVariants[-sel,]}

sel<-which(uniqueVariants$TGP>=0.01)
if(length(sel)>0){uniqueVariants<-uniqueVariants[-sel,]}

print(paste("due", dim(uniqueVariants)[1], sep = "_"))

# FILTER AGAINST SEGDUP

sel<-which(uniqueVariants$SEGDUPS!="")
segdupRem<-uniqueVariants[sel,]
if(length(sel)>0){uniqueVariants<-uniqueVariants[-sel,]}

print(paste("tre", dim(uniqueVariants)[1], sep = "_"))

#sel<-which(uniqueVariants$DBSNP_ID!="" & uniqueVariants$SEGDUPS!="")
#segdupRem<-uniqueVariants[sel,]
#uniqueVariants<-uniqueVariants[-sel,]

# FILTER AGAINST SIFT_PRED REMOVE TOLERATED

# sel<-which(uniqueVariants$SIFT_PRED=="T" & (uniqueVariants$POLY_PRED=="B" | uniqueVariants$POLY_PRED=="P") & (uniqueVariants$MUTASSES_PRED=="N" | uniqueVariants$MUTASSES_PRED=="L"))
# funcRem<-uniqueVariants[sel,]
# uniqueVariants<-uniqueVariants[-sel,]

# sel<-which(uniqueVariants$DBSNP_ID!="" & uniqueVariants$SIFT_PRED=="T" & uniqueVariants$POLY_PRED=="B" & uniqueVariants$MUTASSES_PRED=="N")
# funcRem<-uniqueVariants[sel,]
# uniqueVariants<-uniqueVariants[-sel,]

# REMOVING PUTATIVE POLYMORPHIC GENE

# METHOD1
# removing high frequency mutated genes rsID associated
# <=2 mutations per-gene both rsID associated
# >2 mutations per-gene and 50% rsID associated  

#uniqueVariants<-read.delim(file=paste0(getwd(), "/", fileName), sep="\t", dec=".", stringsAsFactors=F)

# polymorph<-rep(0, dim(uniqueVariants)[1])
# 
# polyrem<-data.frame()
# 
# for(i in 1:length(levels(as.factor(uniqueVariants$GENE)))){
#   
#   sel<-which(uniqueVariants$GENE==levels(as.factor(uniqueVariants$GENE))[i])
#   if(length(sel)>1){
#     if(length(sel)[1]==2){
#       
#       nsnp=length(which(uniqueVariants[sel, 'DBSNP_ID']!=""))
#       if(nsnp==2) {polymorph[sel]=1}
#       
#     } else {
#       
#       nsnp=length(which(uniqueVariants[sel, 'DBSNP_ID']!=""))
#       if(nsnp/length(sel)>=0.50) {polymorph[sel]=1}
#       
#     }
#   }
# }

# METHOD2
# removing high frequency mutated genes rsID associated
# <=2 mutations per-gene both rsID associated
# >2 mutations per-gene and 95% rsID associated  


# polymorph<-rep(0, dim(uniqueVariants)[1])
# 
# polyrem<-data.frame()
# 
# for(i in 1:length(levels(as.factor(uniqueVariants$GENE)))){
#   
#   sel<-which(uniqueVariants$GENE==levels(as.factor(uniqueVariants$GENE))[i])
#   if(length(sel)>1){
#     if(length(sel)[1]==2){
#       
#       nsnp=which(uniqueVariants[sel, 'DBSNP_ID']!="")
#       if(length(nsnp)>0) {polymorph[nsnp]=1}
#       
#     } else {
#       
#       nsnp=which(uniqueVariants[sel, 'DBSNP_ID']!="")
#       if(length(nsnp)>0) {polymorph[nsnp]=1}
#       
#     }
#   }
# }

# DOWNSTREAM TO METHOD 1 OR 2

# CHECK FOR RELEVANT GENES

# pattern<-paste("^", incGenesList[,1], sep="")
# 
# for(i in 1:length(pattern)){
#   chksel<-grep(pattern[i], uniqueVariants$GENE)
#   if(length(chksel)>0) {
#     polymorph[chksel]=0
#   }
# }
# 
# REMOVE POLYMORPHIC

# sel<-which(polymorph==1)
# 
# if(length(sel)>0){
#   
#   polyrem<-uniqueVariants[sel,]
#   
#   uniqueVariants<-uniqueVariants[-sel,]
#   
# }

# GET RSID INFO FROM DBSNP AND REMOVE VARIANTS WITH MAF>0.01
rsIDs<-levels(as.factor(uniqueVariants$DBSNP_ID))
rsIDs<-rsIDs[rsIDs!=""]

# the number of rsIDs that can be passed in to one request is limited to around 600
# therefore, if there are more than 500 I split them in group of 500 and make multiple queries
if(length(rsIDs)>500){
  
  rsIDs_annotated=NULL
  
  n = length(rsIDs)
  groups = (n%/%500)+1
  
  for(i in 1:groups){
    l1 = (500*(i-1))+1
    if(i<groups){l2 = 500*i}else{l2 = n}
    ids = rsIDs[l1:l2]
    rsIDs_annotated<-rbind(rsIDs_annotated, ncbi_snp_query(ids))
  }
  
}else{rsIDs_annotated<-ncbi_snp_query(rsIDs)}

rsIDs_remove<-rsIDs_annotated[which(!is.na(rsIDs_annotated$MAF) & (rsIDs_annotated$MAF>=0.01)), 'Query']
sel<-which(uniqueVariants$DBSNP_ID %in% rsIDs_remove)
if(length(sel)>0){uniqueVariants<-uniqueVariants[-sel,]}

print(paste("quattro", dim(uniqueVariants)[1], sep = "_"))

if(length(rsIDs)>0) {write.table(rsIDs, file="rsIDs.tsv", sep="\t", row.names=F, col.names=F, dec=".", quote=F)}
if(dim(segdupRem)[1]>0) {write.table(segdupRem, file="segdup_removed.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)}
#if(dim(polyrem)[1]>0) {write.table(polyrem, file="polymorphic_removed.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)}
#if(dim(funcRem)[1]>0) {write.table(funcRem, file="function_pred_removed.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)}
if(dim(uniqueVariants)[1]>0) {write.table(uniqueVariants, file="somatic_variants_filtered.tsv", row.names=F, col.names=T, sep="\t", quote=F, dec=".")}

###

# POISSON REGRESSION

# model<-glm(Freq~Var1, data=tab2, family=poisson(link = "log"))
# model<-glm(Freq~Var1, data=tab, family=poisson(link = "identity"))
# pvals<-coef(summary(model))[,4]
# pvals<-pvals[2:length(pvals)]
# pvals<-p.adjust(pvals, method="fdr", n=length(pvals))
# gene_names=substring(names(pvals), 5, length(names(pvals)))
# sigs_gene<-data.frame(gene_names, pvals)
# sel <- which(sigs_gene$pvals<0.05)
# genelist <- sigs_gene[sel, 'gene_names']

###

# tab6<-as.data.frame(table(uniqueVariants$GENE))
# write.table(tab6, file="tab6.tsv", sep="\t", dec=".", row.names=F, col.names=T, quote=F)
# length(intersect(tab6$Var1, tab$Var1))

