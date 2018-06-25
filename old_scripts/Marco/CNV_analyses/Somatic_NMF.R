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

### SETTING WORKING DIRECTORY

#wd = "/media/marco/DATA2/NGS/WES/Projects/Results_AE/signatures"
wd = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/signatures"
od = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/signatures/170519_new_analysis"
setwd(wd)

### LOADING LIBRARY

source("/media/DATA1/DATA/Software/script-R/Functions.R")

library(SomaticSignatures)

library(BSgenome.Hsapiens.UCSC.hg19)

library(survival)

library(RColorBrewer)

library("gplots")


### PREPROCESS VARIANTS LIST TO OBTAIN A VCF FILE

#FILENAME="somatic_variants_filtered_v2.tsv"
FILENAME="170427_snv_signature.tsv"

dset<-read.delim(paste(getwd(), FILENAME, sep="/"), header=T, sep="\t", dec=".")

#dset[which(is.na(dset$X1000g.freq)), "X1000g.freq"]<-0  #Metto a 0 eventuali NA sul campo MAF

sel<-which(dset$CHR=="23")
if(length(sel)>0) dset[sel, 'CHR']="X"

sel<-which(dset$CHR=="24")
if(length(sel)>0) dset[sel, 'CHR']="Y"

dset[, 'CHR']<-paste("chr", dset[, 'CHR'], sep="")

dset$SAMPLE<-paste("Sample_", dset[, 'SAMPLE'], sep="")

#FILTER DATA BY MAF<0.05
#dset<-dset[which(dset$DBSNP.FREQ<0.05 & dset$X1000G.FREQ<0.05),]
#dset<-dset[which(dset$X1000g2012apr_ALL<0.05),]

# SET "chr" FIELD 
#sel<-which((dset$chr!="MT"))
#dset[sel, 'chr']<-paste("chr", dset[sel, 'chr'], sep="")

# REMOVE MT GENES

#sel<-which((dset$CHR=="MT"))
#dset<-dset[-sel, ]


#FILTER BY MUT TYPE !="SYNONYMOUS"
#dset<-dset[which(dset$type.1=="nonsynonymous"), ]

### CONSTRUCTION OF THE VRANGE OBJECT

ir<-IRanges(start=dset$START, end=dset$END)

vr<-VRanges(seqnames=dset$CHR, ranges=ir, ref=dset$REF, alt=dset$ALT, sampleNames=dset$SAMPLE, gene=dset$GENE)


### MUTATION CONTEXT RETRIVAL

mutContext=mutationContext(vr, BSgenome.Hsapiens.UCSC.hg19, k = 3, unify = T)

mutMatrix<-motifMatrix(mutContext, group = "sampleNames", normalize = T)

head(mutMatrix)


### NMF DECOMPOSITION

# Determine the number of sigs

n_sigs = 2

gof_nmf= assessNumberSignatures(mutMatrix, n_sigs, nReplicates = 10)

plotNumberSignatures(gof_nmf)

# Decomposition

sigs_nmf = identifySignatures(mutMatrix, n_sigs, nmfDecomposition)

write.table(sigs_nmf@samples, file=paste0(od, "/", "contrib_2.tsv"), sep="\t", dec=".", quote=F, row.names=T, col.names=NA) 

# Plots

#jpeg("Mutation_spectrum.jpeg", quality = 100, bg = "white", res = 600, width = 20, height = 20, units = "cm")

#plotMutationSpectrum(mutContext, group = "sampleNames")

#dev.off()

jpeg(paste0(od,"/Signatures_2.jpeg"), quality = 100, bg = "white", res = 600, width = 20, height = 20, units = "cm")

plotSignatures(sigs_nmf)

dev.off()

#jpeg("SamplesMap.jpeg", quality = 100, bg = "white", res = 600, width = 16, height = 16, units = "cm")

#plotSampleMap(sigs_nmf)

#dev.off()

jpeg(paste0(od,"/SamplesContrib_2.jpeg"), quality = 100, bg = "white", res = 600, width = 30, height = 30, units = "cm")

plotSamples(sigs_nmf)

dev.off()

### HIERARCHICAL CLUSTERING

# Ward Hierarchical Clustering

d <- dist(sigs_nmf@samples, method = "euclidean") # distance matrix

fit <- hclust(d, method="ward.D")

jpeg(paste0(od,"/SamplesDendro_2.jpeg"), quality = 100, bg = "white", res = 600, width = 16, height = 16, units = "cm")

plot(fit, main="Samples clustered by signatures", xlab="Samples", ylab="Distance", sub="", cex.main=1, cex.lab=0.8, cex.axis=0.8, cex=0.25)

groups <- cutree(fit, k=2) # cut tree into k clusters

# draw dendogram with red borders around the k clusters
rect.hclust(fit, k=2, border="darkblue") 

dev.off()

samples<-unlist(lapply(strsplit(names(groups), "_"), '[[', 2))
groups<-data.frame(samples, groups)
colnames(groups)<-c("SAMPLE", "GROUP")
write.table(groups, file=paste0(od, "/", "groups_2.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)

# STATS





# Plot heatmap

# library(ggplot2)
# 
# base_size <- 9
# 
# (p <- ggplot(as.data.frame(sigs_nmf@samples), aes(variable, Name)) + geom_tile(colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))
 
# p + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
#   + opts(legend.position = "none", axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

#------------------------------------ Heatmap.2

mypalette<-brewer.pal(9, "BuPu")

dd <- function(x) dist(x, method="euclidean")
chclust <- function(d) hclust(d, method="ward.D")

jpeg("SamplesClust.jpeg", quality = 100, bg = "white", res = 600, width = 8, height = 8, units = "cm")

heatmap.2(sigs_nmf@samples, distfun=dd, hclustfun=chclust, Colv=NA, dendrogram="row", key=F, trace="none", cexCol=0.8, cexRow=0.8,
          colsep=F, rowsep=T, sepwidth=c(0.05, 0.05), sepcolor="white") #, col=mypalette
dev.off()


### OVERALL SURVIVAL

SURVFILENAME="/survival.csv"

dset_surv<-read.csv2(paste(getwd(), SURVFILENAME, sep=""), sep = "\t", as.is=T)

STATUS<-rep(0, dim(dset_surv)[1])
STATUS[which(dset_surv$DEATH_DATE!="")]<-1

dset_surv<-cbind(dset_surv, STATUS)

MAX<-max(dset_surv$SURV_MONTHES)

#Between groups
KMfit<-survfit(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, dset_surv)
plot(KMfit, xlab="Months", ylab="Survival", xlim=c(0, 120), bty="n", xaxt="n", main="Between Groups Survival", conf.int=F, mark.time=T, cex.axis=0.8, col=c("black", "",4), lwd=2)
axis(side=1, at = c(0,12, 24, 36, 48, 60, 72, 84, 96, 108, 120), cex.axis=0.8, lwd=2)
axis(side=2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=0.8, lwd=2)
KMfit_log_rank_test<-survdiff(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, dset_surv)
p.val<- 1-pchisq(KMfit_log_rank_test$chisq, length(KMfit_log_rank_test$n) - 1)
legend(MAX-30, 0.9, c("1", "2", "3"), col=c(2,3,4), lty=1, lwd=2, cex=0.8, bty="n")
text(MAX/2, 0.95, paste("p= ", round(p.val, 3), sep=""), cex=0.8)

#Between groups 1-2
surv_set_1_vs_2<-dset_surv[which(dset_surv$GROUP!=3),]
KMfit<-survfit(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_1_vs_2)
plot(KMfit, xlab="Months", ylab="Survival", xlim=c(0, MAX), bty="n", xaxt="n", main="Between Groups Survival", conf.int=F, mark.time=T, cex.axis=0.8, col=c(2,3,4), lwd=2)
axis(side=1, at = c(0,12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132,144), cex.axis=0.8, lwd=2)
axis(side=2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=0.8, lwd=2)
KMfit_log_rank_test<-survdiff(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_1_vs_2)
p.val<- 1-pchisq(KMfit_log_rank_test$chisq, length(KMfit_log_rank_test$n) - 1)
legend(MAX-30, 0.9, c("1", "2"), col=c(2,3,4), lty=1, lwd=2, cex=0.8, bty="n")
text(MAX/2, 0.95, paste("p= ", round(p.val, 3), sep=""), cex=0.8)

#Between groups 1-3
surv_set_1_vs_3<-dset_surv[which(dset_surv$GROUP!=2),]
KMfit<-survfit(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_1_vs_3)
plot(KMfit, xlab="Months", ylab="Survival", xlim=c(0, MAX), bty="n", xaxt="n", main="Between Groups Survival", conf.int=F, mark.time=T, cex.axis=0.8, col=c(2,3,4), lwd=2)
axis(side=1, at = c(0,12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132,144), cex.axis=0.8, lwd=2)
axis(side=2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=0.8, lwd=2)
KMfit_log_rank_test<-survdiff(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_1_vs_3)
p.val<- 1-pchisq(KMfit_log_rank_test$chisq, length(KMfit_log_rank_test$n) - 1)
legend(MAX-30, 0.9, c("1", "3"), col=c(2,3,4), lty=1, lwd=2, cex=0.8, bty="n")
text(MAX/2, 0.95, paste("p= ", round(p.val, 3), sep=""), cex=0.8)

#Between groups 2-3
surv_set_2_vs_3<-dset_surv[which(dset_surv$GROUP!=1),]
KMfit<-survfit(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_2_vs_3)
plot(KMfit, xlab="Months", ylab="Survival", xlim=c(0, MAX), bty="n", xaxt="n", main="Between Groups Survival", conf.int=F, mark.time=T, cex.axis=0.8, col=c(2,3,4), lwd=2)
axis(side=1, at = c(0,12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132,144), cex.axis=0.8, lwd=2)
axis(side=2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=0.8, lwd=2)
KMfit_log_rank_test<-survdiff(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_2_vs_3)
p.val<- 1-pchisq(KMfit_log_rank_test$chisq, length(KMfit_log_rank_test$n) - 1)
legend(MAX-30, 0.9, c("2", "3"), col=c(2,3,4), lty=1, lwd=2, cex=0.8, bty="n")
text(MAX/2, 0.95, paste("p= ", round(p.val, 3), sep=""), cex=0.8)

#Between groups 2-others
surv_set_2<-dset_surv
surv_set_2[which(surv_set_2$GROUP!=2), "GROUP"]<-1
KMfit<-survfit(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_2)

jpeg("OSClust.jpeg", quality = 100, bg = "white", res = 300, width = 10, height = 10, units = "cm")

plot(KMfit, xlab="Months", ylab="Survival", xlim=c(0, 120), bty="n", xaxt="n", main="Overall Survival by signature", conf.int=F, mark.time=T, cex.axis=0.6, cex.main=0.8, cex.lab=0.6, col=c(2,3,4), lwd=2)
axis(side=1, at = c(0,12, 24, 36, 48, 60, 72, 84, 96, 108, 120), cex.axis=0.6, lwd=2)
axis(side=2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=0.6, lwd=2)
KMfit_log_rank_test<-survdiff(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_2)
p.val<- 1-pchisq(KMfit_log_rank_test$chisq, length(KMfit_log_rank_test$n) - 1)
legend(90, 0.9, c("Others", "S3"), col=c(2,3,4), lty=1, lwd=2, cex=0.6, bty="n")
#text(MAX/2, 0.95, paste("p= ", round(p.val, 3), sep=""), cex=0.6)

dev.off()


# At 72 monthes

sel<-which(dset_surv$LAST_FOLLOW_UP_DATE=="")
dset_surv[sel, "LAST_FOLLOW_UP_DATE"]<-dset_surv[sel, "DEATH_DATE"]

at_monthes=72
SURV_MONTHES_fixed<-calc_follow_up_monthes(dset_surv$DIAGNOSIS_DATE, dset_surv$LAST_FOLLOW_UP_DATE, at_monthes)

STATUS_fixed<-rep(0, dim(dset_surv)[1])
STATUS_fixed<-calc_dead_at_monthes(dset_surv$DIAGNOSIS_DATE, dset_surv$DEATH_DATE, at_monthes)
dset_surv<-cbind(dset_surv, SURV_MONTHES_fixed, STATUS_fixed)

#Between groups 2-others
surv_set_2<-dset_surv
surv_set_2[which(surv_set_2$GROUP!=2), "GROUP"]<-1
KMfit<-survfit(Surv(SURV_MONTHES_fixed, STATUS_fixed, type="right")~GROUP, surv_set_2)

jpeg("OSClust.jpeg", quality = 100, bg = "white", res = 300, width = 10, height = 10, units = "cm")

plot(KMfit, xlab="Months", ylab="Survival", xlim=c(0, 120), bty="n", xaxt="n", main="Overall Survival by signature", conf.int=F, mark.time=T, cex.axis=0.6, cex.main=0.8, cex.lab=0.6, col=c(2,3,4), lwd=2)
axis(side=1, at = c(0,12, 24, 36, 48, 60, 72, 84, 96, 108, 120), cex.axis=0.6, lwd=2)
axis(side=2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=0.6, lwd=2)
KMfit_log_rank_test<-survdiff(Surv(SURV_MONTHES, STATUS, type="right")~GROUP, surv_set_2)
p.val<- 1-pchisq(KMfit_log_rank_test$chisq, length(KMfit_log_rank_test$n) - 1)
legend(90, 0.9, c("Others", "S3"), col=c(2,3,4), lty=1, lwd=2, cex=0.6, bty="n")
#text(MAX/2, 0.95, paste("p= ", round(p.val, 3), sep=""), cex=0.6)

dev.off()

### ASSOCIATION ANALYSIS

### ANALYSIS OF PATHWAYS



### ANALYSIS OF MUTATIONS FOR DETECTION OF LOSS OF FUNCTION

### END