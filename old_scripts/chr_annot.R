
########### TAKE A LIST OF GENE SYMBOLS ("symbol" as column name) AND ANNOTATE CHROMOSOMES ############

library(org.Hs.eg.db)

# load gene list
genes = read.delim("/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/170329_chr_toBeAnnotated.txt", header=T, sep="\t", stringsAsFactors = F, check.names = F, na.strings = "na")

# convert to ENTREZ ids
#genes = select(org.Hs.eg.db, keys=genes$symbol, column="ENTREZID", keytype = "SYMBOL")
genes = select(org.Hs.eg.db, keys=genes$symbol, column="ENTREZID", keytype = "ALIAS") # ALIAS retrieves ENTREZ ids in case of 1:many matches

# get genes whose entrez id was not found (NAs)
na = genes[!complete.cases(genes),]
na$chr = rep(NA, dim(na)[1])

# get genes whose entrez id was found
genes = genes[complete.cases(genes),]

# annotate CHR from ENTREZ ids
genes$chr = as.vector(unlist(mget(genes$ENTREZID,envir=org.Hs.egCHR, ifnotfound=NA)))

# re-attach NAs below
genes = rbind(genes, na)

# get list of genes for which multiple ENTREZ ids were found
double = genes[which(duplicated(genes$ALIAS)), "ALIAS"]

# add column for "double" ids
genes$double = rep("n", dim(genes)[1])

# note which genes have multiple alias
for(r in 1:dim(genes)[1]){if(genes$ALIAS[r] %in% double){genes$double[r] = "y"}}

colnames(genes)[1] = "symbol"


write.table(genes, file="/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/170329_chr_annotated.tsv", sep="\t", row.names=F, col.names=T, dec=".", quote=F)
