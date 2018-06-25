#------------------------------------------------------------------------------------------------------------------------
# DEPENDENCIES
#------------------------------------------------------------------------------------------------------------------------
source("/Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/functions.R")  


library(org.Hs.eg.db)
library(annotate)
library(ggplot2)

amltp53DIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/tp53"


#------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------- FUNCTIONS --------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------
#-------------------# 170404 - GENES ON MULTIPLE CHROMOSOMES (ACCORDING TO NEXUS OUTPUT) #--------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

### Determine which genes in our NEXUS output are indicated to lie in multiple chromosomes

# load and edit nexus output
nexus = read.delim(paste0(amltp53DIR, "/nexus.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
colnames(nexus)[length(colnames(nexus))] = "symbol"

# extract CHR info from column with genomic coordinates
nexus$chr = sapply(nexus[,2], function(x) strsplit(x, ":")[[1]][1])

# remove LOH events, allelic imbalance, non-genic events
nexus = nexus[nexus$Event!="LOH",]
nexus = nexus[nexus$Event!="Allelic Imbalance",]
nexus = nexus[nexus$symbol!="",]

# extract each single gene and make a table with sample and chromosomes info
cgs=NULL
for(r in 1:dim(nexus)[1]){
  if(r%%100==0){cat(r,"/",dim(nexus)[1],"\n")}
  chr = nexus$chr[r]
  genes = unlist(strsplit(nexus$symbol[r], ", "))
  patient = nexus$Sample[r]
  cgs = rbind(cgs, data.frame(sample=rep(patient, length(genes)), chr=rep(chr, length(genes)), gene=genes, stringsAsFactors=F))}

# remove duplicate rows from 3 columns table
cgs = unique(cgs)
# get only "gene" and "chr" columns
cg = unique(cgs[,c(2,3)])
# sort by "gene"
cg = cg[order(cg$gene),]
# select genes that are duplicated in this table (which means they are associated to more than one chromosome)
sel = unique(cg[which(duplicated(cg$gene)), "gene"])
# subset 3 column table for genes that appear in the duplicated genes selection
cgs_sel = cgs[cgs$gene %in% sel,]

write.table(sel, file=paste0(amltp53DIR, "/170404_genes_on_multiple_chromosomes.txt"), sep="\t", dec=".", col.names=F, row.names=F, quote=F)
write.table(cgs_sel, file=paste0(amltp53DIR, "/170404_genes_on_multiple_chromosomes_with_pat_and_chr.tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)



#---------------------------------------------------------------------------------------------------------------------
#-------------------------------------------# 2016 - NEXUS SETTING_1 [WT VS. MUT] #--------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/1611_tp53/setting1_wt_mut/")

# transform DOS based file to UNIX based (Mac) and edit in order to be compatible with R
system("sh /Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/1611_tp53_editNexus.sh < classic.txt > classic_ed.txt")
system("dos2unix -c Mac classic_annot.txt")


# import data
classic = read.table("classic_ed.txt", header = T, sep="\t")
annot = read.table("classic_annot.txt", header = T, sep="\t")
merged = unique(merge(classic,annot, by="Gene.Symbol"))

write.table(merged, file="merged.txt", sep="\t", row.names=FALSE, quote=FALSE)

# filter for significants after multiple testing
merged1 = merged[merged$q_bound < 0.05,]

# remove double genes
to_delete = vector()
for(i in 2:dim(merged1)[1]){
  if(merged1$Gene.Symbol[i] == merged1$Gene.Symbol[i-1]){
    to_delete = append(to_delete,i)
  }
}
merged2 = merged1[-to_delete,]

#remove antisense, linc, miRNA, orf, LOC
merged3 = merged2[grepl("LINC",merged2$Gene.Symbol) == FALSE,]
merged3 = merged3[merged3$miRNAs == ".",]
merged3 = merged3[grepl("orf",merged3$Gene.Symbol) == FALSE,]
merged3 = merged3[grepl("-AS",merged3$Gene.Symbol) == FALSE,]
merged3 = merged3[grepl("MIR",merged3$Gene.Symbol) == FALSE,]
merged3 = droplevels(merged3[grepl("LOC",merged3$Gene.Symbol) == FALSE,])

write.table(merged3, file="merged_filtered.txt", sep="\t", row.names=FALSE, quote=FALSE)


#####################################
# pathway analysis with Marco's script
#####################################

# extract gene symbols list and write to file
gene_symbol_list = unique(merged3$Gene.Symbol)
write("GENE", file="genes_list.txt")
for(i in 1:length(gene_symbol_list)){
  write(as.character(gene_symbol_list[i]), file="genes_list.txt", append=T)
}

# run dbMapping script, set new WD and edit reactome output list
system("Rscript /Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/DBmapping.R < genes_list.txt")
setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/1611_tp53/setting1_wt_mut/db_mapping_output/")
system("sh /Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/1611_tp53_editReactome.sh < reactome_list.txt > reactome_list_ed.txt")

# import edited reactome list, add "counts" and "Gene.Symbol" columns
# "counts" expresses how many genes share each pathway
reactome = read.table("reactome_list_ed.txt", header = T, sep="\t", quote="")
counts = factor(1:dim(reactome)[1])
Gene.Symbol = factor(1:dim(reactome)[1])
reactome$counts = counts
reactome$counts = 0
reactome$Gene.Symbol = "."
gene_pathname = unique(reactome[,c(1,2)])
counts_table = table(gene_pathname$PATHNAME)
for(i in 1:dim(reactome)[1]){
  if(i%%10 == 0){
    message(i/dim(reactome)[1])
  }
  reactome$Gene.Symbol[i] = lookUp(as.character(reactome$ENTREZID[i]), 'org.Hs.eg', 'SYMBOL')
  if(reactome$PATHNAME[i] == "."){
    reactome$counts[i] = "."
  } else{
    if(reactome$PATHNAME[i] %in% names(counts_table)){
      symbol = reactome$PATHNAME[i]
      reactome$counts[i] == as.integer(counts_table[symbol])
    }
  }
}


######
rimerged = merge(merged,reactome, by="Gene.Symbol")
#######

write.table(reactome, file="reactome_path_counts.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(rimerged, file="rimerged.txt", quote=F, row.names=F, col.names=T, sep="\t")






#----------------------------------------------------------------------------------------------------------------
#-------------------------------------------# 2016 - NEXUS SETTING_2 #--------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

MAINPATH = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/1611_tp53/setting2_wt_mut_mutDel"
SCRIPTPATH = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/"

MAINPATH = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/1611_tp53/setting2_prova2"
SCRIPTPATH = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/"
comparisons = c("mutDel_wt", "mutMutDel_wt", "mut_wt", "mut_mutDel")

mergeFilter <- function(folder, MAINPATH, SCRIPTPATH){ 
  
  # set paths and directories
  SCRIPTPATH = SCRIPTPATH
  FOLDERPATH = file.path(MAINPATH, folder)
  DBPATH = file.path(FOLDERPATH, "db_mapping_output")
  prefix = paste0("124_", folder, "_")
  dir.create(DBPATH)
  
  setwd(FOLDERPATH)

  # transform DOS based file to UNIX based (Mac) and edit in order to be compatible with R
  files = list.files(FOLDERPATH)
  for(i in 1:length(files)){
    system(paste0("dos2unix -c Mac ", files[i]))}
  
  system(paste0("sh ", SCRIPTPATH, "editNexus.sh ", FOLDERPATH))

  # import data
  types = (c("classic","combined","peaks"))
  for(i in 1:length(types)){
    print(types[i])
    intervals = read.table(paste0(prefix, types[i], "_ed.txt"), header = T, sep="\t", quote = "", comment.char = "")
    annot = read.table(paste0(prefix, types[i], "_annot_ed.txt"), header = T, sep="\t", quote = "", comment.char = "")
    merged = unique(merge(intervals,annot, by="Gene.Symbol"))
    
    write.table(merged, file=paste0(types[i], "_merged.txt"), sep="\t", row.names=FALSE, quote=FALSE)

    # filter for significants after multiple testing
    merged1 = droplevels(merged[merged[,13] < 0.05,])
    # remove double genes
    if(!(dim(merged1)[1] < 1)){
      to_delete = vector()
      for(j in 2:dim(merged1)[1]){
        if(merged1[,1][j] == merged1[,1][j-1]){
          to_delete = append(to_delete,j)
        }
      }
      print(dim(merged1))
      print(is.null(length(to_delete)))
      if(!(is.null(length(to_delete)))){
        merged1 = droplevels(merged1[-to_delete,])
        print(dim(merged1))}
      
      write.table(merged1, file=paste0(types[i], "_merged_filtered.txt"), sep="\t", row.names=FALSE, quote=FALSE)
      
      #####################################
      # pathway analysis with Marco's script
      #####################################
      
      # extract gene symbols list and write to file
      gene_symbol_list = unique(merged1$Gene.Symbol)
      write("GENE", file=paste0(DBPATH, "/", types[i], "_genes_list.txt"))
      for(j in 1:length(gene_symbol_list)){
        write(as.character(gene_symbol_list[j]), file=paste0(DBPATH, "/", types[i], "_genes_list.txt"), append=T)
      }
      
      # run dbMapping script and edit reactome output list
      system(paste0("Rscript ", SCRIPTPATH, "DBmapping.R ", DBPATH, "/", types[i], "_genes_list.txt"))
      system(paste0("sh ", SCRIPTPATH, "editReactome.sh < ", DBPATH, "/REACTOMEList.", types[i], "_genes_list.tsv > ", DBPATH, "/", types[i], "_reactome_ed.txt"))
      
      # import edited reactome list, add "counts" and "Gene.Symbol" columns
      # "counts" expresses how many genes share each pathway
      reactome = read.table(paste0(DBPATH, "/", types[i], "_reactome_ed.txt"), header = T, sep="\t", quote="")
      counts = factor(1:dim(reactome)[1])
      Gene.Symbol = factor(1:dim(reactome)[1])
      reactome$counts = counts
      reactome$counts = 0
      reactome$Gene.Symbol = "."
      gene_pathname = unique(reactome[,c(1,2)])
      counts_table = table(gene_pathname$PATHNAME)
      for(j in 1:dim(reactome)[1]){
        if(j%%100 == 0){
          message(j/dim(reactome)[1])}
        if(is.na(reactome$ENTREZID[j])){
          reactome$ENTREZID[j] = "."}
        reactome$Gene.Symbol[j] = lookUp(as.character(reactome$ENTREZID[j]), 'org.Hs.eg', 'SYMBOL')
        if(reactome$PATHNAME[j] == "."){
          reactome$counts[j] = "."
        } else {
          if(reactome$PATHNAME[j] %in% names(counts_table)){
            symbol = reactome$PATHNAME[j]
            reactome$counts[j] = as.integer(counts_table[symbol])}
        }
      }
      reactome[,5] = unlist(reactome[,5])
      
      write.table(reactome, file=paste0(DBPATH, "/", types[i], "_reactome_path_counts.txt"), quote=F, row.names=F, col.names=T, sep="\t")
    }
  }
}

for(i in 1:length(comparisons)){
  mergeFilter(comparisons[i], MAINPATH, SCRIPTPATH)
}




######
rimerged = merge(merged,reactome, by="Gene.Symbol")
write.table(rimerged, file="rimerged.txt", quote=F, row.names=F, col.names=T, sep="\t")
#######




#-----------------------------------------------------------------------------------------------------------------
#------------------------------------# 2016 - PATHWAY ENRICHMENT - REACTOME #--------------------------------------------
#----------------------------------------------------------------------------------------------------------------

setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/1611_tp53/setting2_prova2/mut_wt/db_mapping_output")
biocLite("ReactomePA")
library(ReactomePA)
df = read.table("classic_genes_list.txt", header = T, sep="\t")
geneSymb = as.character(df[,1])
a = org.Hs.egSYMBOL2EG
symbol2entrez = as.list(a[mappedkeys(a)])

entrezID = character()
for(i in 1:length(geneSymb)){
  entrez = symbol2entrez[[geneSymb[i]]]
  if(is.null(entrez) == TRUE){
    print(geneSymb[i])
  }else{
    entrezID = append(entrezID, symbol2entrez[[geneSymb[i]]])
  }
}

pa = enrichPathway(gene=entrezID, organism="human", pvalueCutoff=0.05, readable=T)
head(as.data.frame(pa))

barplot(pa)

#####################################
#####################################
#####################################


# filter for significants after multiple testing
classic1 = classic[classic$q_bound < 0.05,]

classic2 = droplevels(classic1[classic1$Gene.Symbols != ".",])

genes_cl = sort(unique(classic2$Gene.Symbols))
genes_annot = sort(unique(annot$Gene.Symbol))


write(levels(genes_cl), file="ALL/1611_tp53/enriched_genes.txt", ncolumns=length(genes_cl), sep=" ")




#####################################
name = factor(c("gino","pina","lino","rino","rino","john"))
age = factor(c(23,45,2,97,97,"."))
hair = factor(c("b","b","b","w","w","w"))
name_age = data.frame(name,age,hair)
name = factor(c("gino","pina","lino","rino","john"))
sex = factor(c("m","f","m","m","m"))
hair = factor(c("b","b","b","w","w"))
name_sex = data.frame(name,sex,hair)
p = merge(name_age,name_sex)
unique(merge(name_age,name_sex))
#####################################

for(i in 2:dim(annot)[1]){
  if(annot$Gene.Symbol[i] == annot$Gene.Symbol[i-1]){
    print(annot$Gene.Symbol[i])
  }
}
#####################################

comparison = as.character(m2$chr.x) == as.character(m2$chr.y)
for(i in 1:length(comparison)){
  if(comparison[i] == "FALSE"){
    print(m2[i,"Gene.Symbol"])
    print(m2[i,"chr.x"])
    print(m2[i,"chr.y"])
    print("------------")
  }
}

#####################################
for(i in 1:length(genes_annot)){
  print(genes_cl[i], genes_annot[i])
  #  if(genes_cl[i] != genes_annot[i])
  #    print(i)
}
#####################################

# to transform NA to "." in a data frame
reactome = as.matrix(reactome)
reactome[is.na(reactome)] = "."
reactome = as.data.frame(reactome)




#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------- FOR TESTING PURPOSE --------------------------------------------
#-------------------------------------------------------------------------------------------------------------

del = read.delim(paste0(amltp53DIR, "/arrayTable_del.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
dup = read.delim(paste0(amltp53DIR, "/arrayTable_dup.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
gain = read.delim(paste0(amltp53DIR, "/arrayTable_gain.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
loss = read.delim(paste0(amltp53DIR, "/arrayTable_loss.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)


r = 11



sym = nexus$symbol
head(sym)

sym[grep("DEC", sym)]
nexus[48689,]

geneTable = del
colnames(geneTable)[grep("SEP", colnames(geneTable))]
colnames(geneTable)[sel]

loss[,"SEPT9"]

nexus[nexus[,"symbol"] %in% date_gene[,"gene"], "symbol"]
nexus2[nexus2[,"symbol"] %in% date_gene[,"date"], "symbol"]

nexus[nexus[,"symbol"] %in% date_gene[,"date"], "symbol"] = sapply(nexus[nexus[,"symbol"] %in% date_gene[,"date"], "symbol"], function(x) x=date_gene[date_gene[,"date"]==x,"gene"])

