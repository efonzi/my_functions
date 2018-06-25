
# This script must be run in a folder containing one or more subfolders whose names must contain 4 strings separated by "_":
# 
# (1) name of the platform, i.e. "array"
# (2) name of the phenotypic character upon which the sample classification is done, i.e. "chromothripsis"
# (3) name of the event, i.e. "dup" (this could be double, like "gain-dup")
# (4) numbers identifying the cases (left) and ctrls (right)
#
# Example: "array_chromothripsis_dup_1vs2"
#
#
# The subfolders must contain a table with names of samples and numbers indentifying cases and ctrls.
# The name of this file must contain the string ‘arrayNames'. There must not be other files that contain this string in their name.
#
# Example:
#
# Sample tp53  chromothripsis  <--- table titles
# sample1   1     2
# sample2   1     1
# sample3   2     2
# sample4   2     1
# sample5   1     1
# sample6   2     2
#
#
#
# The subfolders must also contain one or two PATIENTxGENE matrices with the counts of a specific event.
# The...............
#


#----------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------- DEPENDENCIES ----------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

tryCatch({source("/Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/functions.R") 
}, error = function(err){
  print("can't find --> /Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/functions.R")
  print("trying with --> /home/PERSONALE/eugenio.fonzi2/script/functions.R")
}, finally = {source("/home/PERSONALE/eugenio.fonzi2/script/functions.R")})

tryCatch({source("/Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/functions_gene_pathway_enrichment.R") 
}, error = function(err){
  print("can't find --> /Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts/functions_gene_pathway_enrichment.R")
  print("trying with --> /home/PERSONALE/eugenio.fonzi2/script/functions_gene_pathway_enrichment.R")
}, finally = {source("/home/PERSONALE/eugenio.fonzi2/script/functions_gene_pathway_enrichment.R")})

library(org.Hs.eg.db)
library(reactome.db)
library(clusterProfiler)
library(ReactomePA)
library(GO.db)
library(KEGGREST)

#----------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------- LOOP THROUGH FOLDERS  ----------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
setwd(wd)

#args <- args[-1] #avoid problems if "cutoff" is in the directory name

folders = list.files(getwd())

pathDIR = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/pathway_databases"

cat("\nfolders are\n")
print(folders)
cat("\n\n arguments are\n")
print(args)

#----------------------------- NEXUS to 4 MATRICES -------------------------

if("-n2m" %in% args){
   
    cat("nexus2matrices\n")

    # load and edit sample list
    sampleID = read.delim(list.files(getwd(), 'arrayNames'), header=T, sep="\t", stringsAsFactors = F, check.names = F)
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))

    # load and edit nexus output
    nexus = read.delim("nexus.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F)
    nexus$Sample<-unlist(lapply(strsplit(nexus$Sample, " "), "[[", 1))

    # select desired samples in nexus output
    nexus = nexus[which(nexus$Sample %in% sampleID$Sample),]

    # edit column name
    colnames(nexus)[length(colnames(nexus))] = "symbol"

    # remove LOH events, non-genic events and change other event names
    nexus = nexus[nexus$Event!="LOH",]
    nexus = nexus[nexus$Event!="Allelic Imbalance",]
    nexus = nexus[nexus$symbol!="",]
    nexus[nexus$Event=="CN Loss", "Event"] = "loss"
    nexus[nexus$Event=="CN Gain", "Event"] = "gain"
    nexus[nexus$Event=="High Copy Gain", "Event"] = "dup"
    nexus[nexus$Event=="Homozygous Copy Loss", "Event"] = "del"
    
    # extract chromosome info
    nexus$chr = sapply(nexus$`Chromosome Region`, function(x) strsplit(x, ":")[[1]][1])

    # keep only useful columns
    nexus = nexus[,c(1,3,10,11)]
    
    nexus2=NULL
    for(r in 1:dim(nexus)[1]){
      print(r)
      genes = unlist(strsplit(nexus$symbol[r], ", "))
      n=NULL
      for(g in 1:length(genes)){n = rbind(n, nexus[r,])}
      n$symbol = genes
      nexus2 = rbind(nexus2, n)}
    nexus = nexus2
    rm(nexus2)
    
    ####################################################################################
    ## some gene symbols are converted to date by EXCEL - this is specific for amlTP53
    ####################################################################################
    # create conversion table
    date_gene = data.frame(date = c("1-Dec", "1-Mar", "14-Sep", "3-Mar", "6-Sep", "8-Mar", "9-Sep"), gene = c("DEC1", "MARCH1", "SEPT14", "MARCH3", "SEPT6", "MARCH8", "SEPT9"), stringsAsFactors = F)
    # find in "nexus" genes converted to date and convert them back using the conversion table
    nexus[nexus[,"symbol"] %in% date_gene[,"date"], "symbol"] = sapply(nexus[nexus[,"symbol"] %in% date_gene[,"date"], "symbol"], function(x) x=date_gene[date_gene[,"date"]==x,"gene"])
 
    #######################################################################################
    ## nexus assigns some genes to multiple chromosomes - this step is specific for amlTP53
    #######################################################################################
    # load and edit table with details on what genes to keep and on which are the correct chromosomes
    toKeep = read.delim("170404_genes_on_multiple_chromosomes_to_keep.txt", header=T, sep="\t", stringsAsFactors = F, check.names = F)
    toKeep$chr = rep(".", dim(toKeep)[1])
    toKeep[grep("q",toKeep$cytoband), "chr"] = sapply(toKeep[grep("q",toKeep$cytoband), "cytoband"], function(x) paste0("chr", strsplit(x, "q")[[1]][1]))
    toKeep[grep("p",toKeep$cytoband), "chr"] = sapply(toKeep[grep("p",toKeep$cytoband), "cytoband"], function(x) paste0("chr", strsplit(x, "p")[[1]][1]))
    rownames(toKeep) = toKeep$gene
    # loop over list and remove unwanted genes and genes annotated on wrong chromosomes
    for(r in 1:dim(toKeep)[1]){
      gene = toKeep$gene[r]
      chr = toKeep$chr[r]
      keep = toKeep$keep[r]
      if(keep=="no"){nexus = nexus[nexus$symbol!=gene,] # remove genes marked as "no" in "toKeep"
      }else{nexus = nexus[-(which((nexus$symbol==gene & nexus$chr!=chr))),]}} # for genes marked as "yes", keep only events annotated on the correct chromosomes 
    
    # # get unique list of genes with at least an event
    # symbols=NULL
    # for(r in 1:dim(nexus)[1]){
    #   if(r%%100==0){cat("extracting genes for event", r, "/", dim(nexus)[1], "\n")}
    #   symbols = append(symbols, unlist(strsplit(nexus$symbol[r], ", ")))}
    # symbols = sort(unique(symbols))

    # get unique list of genes
    symbols = sort(unique(nexus$symbol))
    # get unique list of patients
    patients = sort(unique(nexus$Sample))

    # create PATIENTxGENE matrices for "del", "dup", "gain", "loss"
    # patients as rownames and first column
    # columns for genes are filled with zeroes
    matrix = data.frame(Sample=patients, matrix(data=0, nrow=length(patients), ncol=length(symbols)))
    colnames(matrix)[-1] = symbols
    rownames(matrix) = patients
    del = matrix
    dup = matrix
    gain = matrix
    loss = matrix
    
    # loop over events and increase by one relevant fields in the matrices
    for(r in 1:dim(nexus)[1]){
      if(r%%100==0){cat("filling PATIENTxGENE matrix for event", r, "/", dim(nexus)[1], "\n")}
      alteration = nexus$Event[r]
      patient = nexus$Sample[r]
      gene = nexus$symbol[r]
      if(alteration == "del"){del[patient, gene] = del[patient, gene] + 1}
      if(alteration == "dup"){dup[patient, gene] = dup[patient, gene] + 1}
      if(alteration == "gain"){gain[patient, gene] = gain[patient, gene] + 1}
      if(alteration == "loss"){loss[patient, gene] = loss[patient, gene] + 1}
    }

    cat("\n writing matrices to file\n")
    write.table(del, file="arrayTable_del.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
    write.table(dup, file="arrayTable_dup.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
    write.table(gain, file="arrayTable_gain.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)
    write.table(loss, file="arrayTable_loss.tsv", sep="\t", dec=".", col.names=T, row.names=F, quote=F)

}

#----------------------------- FOLDER LOOP - TABLE EDIT, FILTER, MERGE AND CONVERT TO PATIENTxGENE MATRIX -------------------------

if("-t2m" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("table2matrix for", folderName, "\n")
    
    # in case of SNP array data
    if(platform == "array"){
      
      # load and edit samples list (remove anything written after a space --> "copy")
      sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=F, check.names=F)
      sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
      
      # if two events have to be merged
      if(grepl("-", event)){
        event1 = strsplit(event, "-")[[1]][1]
        event2 = strsplit(event, "-")[[1]][2]
        geneTable1 = read.delim(paste0("arrayTable_", event1,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
        geneTable2 = read.delim(paste0("arrayTable_", event2,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
        
        # edit sample names (remove anything written after a space)
        geneTable1$Sample<-unlist(lapply(strsplit(geneTable1$Sample, " "), "[[", 1))
        geneTable2$Sample<-unlist(lapply(strsplit(geneTable2$Sample, " "), "[[", 1))
        
        # import sample list, filter out unwanted genes and patients from CHROMOTHRIPSIS dataset
        if(study == "chromot"){
          geneTable1 = filterOutGenesAndPatientsChromothripsis(geneTable1, sampleID)
          geneTable2 = filterOutGenesAndPatientsChromothripsis(geneTable2, sampleID)}
        
        # merge the two tables
        geneTable = mergePatientGeneTables(geneTable1, geneTable2)
        
        # re-add patients filtered-out because they had no events (all their row will be of 0s)
        add = sampleID[which(!sampleID$Sample %in% geneTable$Sample),1] #check samples list for missing samples in the table
        df = data.frame(Sample = add, matrix(data=0, nrow=length(add), ncol=dim(geneTable)[2]-1 ), stringsAsFactors = F) # create dataframe of zeroes and row names of samples to add
        colnames(df) = colnames(geneTable) # adjust column names of dataframe
        geneTable = rbind(geneTable, df) # merge
        
        # quality check: interrupt if any sample name in samples list is not present in PATIENTxGENE matrix
        stopifnot(!FALSE %in% (sampleID$Sample %in% geneTable$Sample))

        # remove genes (columns) with no events (all zeroes)
        sel = which(apply(geneTable[,2:dim(geneTable)[2]], 2, sum) == 0)+1
        if(length(sel)>0){geneTable = geneTable[,-sel]} #+1 is needed to avoid "Sample" column to be removed
        
        write.table(geneTable, file=paste0("matrixPatientGene_", folderName,".tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)}
      
      # if there is only one event and no merging is needed
      else{
        geneTable = read.delim(paste0("arrayTable_", event,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
        
        # edit sample names (remove anything written after a space)
        geneTable$Sample<-unlist(lapply(strsplit(geneTable$Sample, " "), "[[", 1))
        
        # filter out unwanted genes and patients from CHROMOTHRIPSIS dataset
        if(study == "chromot"){geneTable = filterOutGenesAndPatientsChromothripsis(geneTable, sampleID)}
        
        # re-add patients filtered-out because they had no events (all their row will be of 0s)
        add = sampleID[which(!sampleID$Sample %in% geneTable$Sample),1] #check samples list for missing samples in the table
        df = data.frame(Sample = add, matrix(data=0, nrow=length(add), ncol=dim(geneTable)[2]-1 ), stringsAsFactors = F) # create dataframe of zeroes and row names of samples to add
        colnames(df) = colnames(geneTable) # adjust column names of dataframe
        geneTable = rbind(geneTable, df) # merge
        
        # quality check: interrupt if any sample name in samples list is not present in PATIENTxGENE matrix
        stopifnot(!FALSE %in% (sampleID$Sample %in% geneTable$Sample))
        
        # remove genes (columns) with no events (all zeroes)
        geneTable = geneTable[,-(which(apply(geneTable[,2:dim(geneTable)[2]], 2, sum) == 0)+1)] #+1 is needed to avoid "Sample" column to be removed
        
        write.table(geneTable, file=paste0("matrixPatientGene_", folderName,".tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)}
    }
    
    # in case of WES data
    else{
      
      # load event table and samples list
      wesTable = read.delim(paste0("wesTable_", event, ".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      sampleID = read.delim(paste0("wesNames.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      
      # convert gene "ZNF595,ZNF718" to "ZNF595"
      if(study == "aneupl" | study == "chromot"){wesTable[wesTable$GENE=="ZNF595,ZNF718", 6] = "ZNF595"}
      
      # convert wes table to PATIENTxGENE matrix
      cat("converting wes table to PATIENTxGENE matrix \n")
      geneTable = convertWesTableToPatientGeneTable(wesTable)
      
      # edit sample names (remove anything written after a space)
      geneTable$Sample<-unlist(lapply(strsplit(geneTable$Sample, " "), "[[", 1))
      sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
      
      # quality check: interrupt if any sample name in samples list is not present in PATIENTxGENE matrix
      stopifnot(!FALSE %in% (sampleID$Sample %in% geneTable$Sample))
      
      # remove genes (columns) with no events (all zeroes)
      geneTable = geneTable[,-(which(apply(geneTable[,2:dim(geneTable)[2]], 2, sum) == 0)+1)] #+1 is needed to avoid "Sample" column to be removed
      
      write.table(geneTable, file=paste0("matrixPatientGene_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      
    }
    
  }
  
}

#----------------------------- FOLDER LOOP - ANNA - TABLE EDIT, FILTER, MERGE AND CONVERT TO PATIENTxGENE MATRIX -------------------------

if("-t2mANNA" %in% args){
  
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("\n table2matrix for", folderName, "\n")
    
    # load and edit samples list (remove anything written after a space --> "copy")
    sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=F, check.names=F)
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    # if multiple events have to be merged
    if(grepl("-", event)){
      
      ## for each event, load specific "arrayTable" and edit sample names (remove anything written after a space)
      # get all events
      events = unlist(strsplit(event, "-"))
      
      cat("merging and filtering", events, "\n")
      
      # load table for first event, edit sample names and get patients list
      geneTable = read.delim(paste0("arrayTable_", events[1],".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      geneTable$Sample = sapply(geneTable$Sample, function(x) strsplit(x, " ")[[1]][1])
      patients = geneTable$Sample
      
      # loop over other events and load specific table, edit sample names and merge to the previous ones
      for(e in 2:length(events)){
        cat("merging", events[e], "\n")
        matrix = read.delim(paste0("arrayTable_", events[e],".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
        matrix$Sample = sapply(matrix$Sample, function(x) strsplit(x, " ")[[1]][1])
        patients2 = matrix$Sample
        if(FALSE %in% (patients == patients2)){cat("patients issue\n")}
        stopifnot(!FALSE %in% (patients == patients2))
        geneTable = data.frame(Sample=patients, geneTable[,-1] + matrix[,-1], stringsAsFactors=F, check.names = F)}
      
      # filter out genes on Black List
      geneTable = filterOutGenesAnna(geneTable, sampleID)

      # quality check: interrupt if any sample name in samples list is not present in PATIENTxGENE matrix
      stopifnot(!FALSE %in% (sampleID$Sample %in% geneTable$Sample))
      
      # remove genes (columns) with no events (all zeroes)
      sel = which(apply(geneTable[,2:dim(geneTable)[2]], 2, sum) == 0) +1 #+1 is needed to avoid "Sample" column to be removed
      if(length(sel>0)){geneTable = geneTable[,-sel]}
      
      write.table(geneTable, file=paste0("matrixPatientGene_", folderName,".tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)}
    
    # if there is only one event and no merging is needed
    else{
      
      cat("filtering", event, "\n")
      
      # load table and edit sample names (remove anything written after a space)
      geneTable = read.delim(paste0("arrayTable_", event,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      geneTable$Sample<-unlist(lapply(strsplit(geneTable$Sample, " "), "[[", 1))
      
      # filter out genes on Black List
      geneTable = filterOutGenesAnna(geneTable, sampleID)
      
      # quality check: interrupt if any sample name in samples list is not present in PATIENTxGENE matrix
      stopifnot(!FALSE %in% (sampleID$Sample %in% geneTable$Sample))
      
      # remove genes (columns) with no events (all zeroes)
      geneTable = geneTable[,-(which(apply(geneTable[,2:dim(geneTable)[2]], 2, sum) == 0)+1)] #+1 is needed to avoid "Sample" column to be removed
      
      write.table(geneTable, file=paste0("matrixPatientGene_", folderName,".tsv"), sep="\t", dec=".", col.names=T, row.names=F, quote=F)}
    
  }
  
}

#------------------------------------------------ FOLDER LOOP - TRANSFORM MULTIPLE EVENTS TO 1 IN PATIENTxGENE MATRIX ------------------------------------------------

if("-to1" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix\n")
    geneTable = read.delim(paste0("matrixPatientGene_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # change to 1 all values higher than 1
    geneTable = data.frame(Sample=geneTable[,1], apply(geneTable[-1], 2, function(x) sapply(x, function(y) if(y>1){y=1}else{y=y})), stringsAsFactors = F, check.names = F)
    
    write.table(geneTable, file=paste0("matrixPatientGene_allOnes_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  
  }
}
    
#------------------------------------------------ FOLDER LOOP - FISHER ON GENES ------------------------------------------------

if("-gf" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("\nFISHER ON GENES - ", folderName, "\n")
    
    # load and edit relevant samples list (remove anything written after a space --> "copy")
    cat("loading and editing sample table\n")
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
    
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    # create list of cases, ctrls and NAs according to the specific classification criterion
    selCase = separateCaseCtrl(criterion, case, ctrl, sampleID)[[1]]
    selCtrl = separateCaseCtrl(criterion, case, ctrl, sampleID)[[2]]
    selNA = separateCaseCtrl(criterion, case, ctrl, sampleID)[[3]]
    
    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix\n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # patient number checkpoint
    cat(dim(selCase)[1], "cases", dim(selCtrl)[1], "ctrls", dim(selNA)[1], "NAs --->", dim(selCase)[1]+dim(selCtrl)[1]+dim(selNA)[1], "total \n")
    cat(dim(geneTable)[1], "samples in event matrix \n")
    
    ###########################

    Rmatrix = geneCaseCtrlFisher(geneTable, selCase, selCtrl)
    write.table(Rmatrix, file=paste0("geneFisher_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
  
}


#--------------------------------------- FOLDER LOOP - PATHWAY ENRICHMENT (OVER-REPRESENTATION TEST) -----------------------------------------------

if("-or" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("\nPATHWAY ENRICHMENT ANALYSIS - OVER-REPRESENTATION TEST - ", folderName, "\n")
    
    # load and edit relevant samples list (remove anything written after a space --> "copy")
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
    
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix \n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # perform over-representation hypergeometic test on pathways for all patients and concatenate all results
    
    allDatabases = hypergeometricTestOnMultipleIndividuals(geneTable)
    
    reactomeALL = allDatabases[[1]]
    
    keggALL = allDatabases[[2]]
    
    gobpALL = allDatabases[[3]]
    
    write.table(reactomeALL, file=paste0("over-representationREACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(keggALL, file=paste0("over-representationKEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(gobpALL, file=paste0("over-representationGOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
  
}



#######################################################################
#######################################################################
##### TEMPORARY FOR MCF!!!!!!!!!!!
#######################################################################
#--------------------------------------- FOLDER LOOP - PATHWAY ENRICHMENT (OVER-REPRESENTATION TEST) -----------------------------------------------
hypergeometricTestOnMultipleIndividualsMCF <- function(geneTable){
  
  reactomeALL=NULL
  
  
  # loop over table rows (patients)
  for(r in 1:dim(geneTable)[1]){
    
    # get name of the patient
    patient = geneTable[r,1]
    
    # select altered genes for the patient
    genes = colnames(geneTable)[which(geneTable[r,] == 1)]
    
    cat("     \n")
    cat("PATIENT", r, "of", dim(geneTable)[1], " - ", length(genes), "events", " - ", folderName, "\n")
    
    if(length(genes)>0){
      
      cat("pathway over-representation test \n")
      
      # transform gene symbols to entrez ids
      ids = as.vector(unlist(mget(genes,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
      ids = ids[complete.cases(ids)]
      
      # perform over-representation test for reactome, add sample name, add case/ctrl classification number and
      # concatenate to "reactomeALL"
      if(length(ids)>0){
        if(TRUE %in% (ids %in% keys(reactome.db, keytype="ENTREZID"))){ # to check whether the ids map in reactome
          reactome = try(enrichPathway(ids, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1))
          if(!inherits(reactome, "NULL")){
            reactome = reactome@result
            reactome$Sample = patient
            reactomeALL = rbind(reactomeALL, reactome)
            cat("reactome --> DONE \n")
          }else{cat("SKIP reactome \n")}
        }else{cat("SKIP reactome \n")}}
      
      
    }
  }
  
  return(reactomeALL)
  
}

if("-orMCF" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("\nPATHWAY ENRICHMENT ANALYSIS - OVER-REPRESENTATION TEST - ", folderName, "\n")

    # load and edit relevant samples list (remove anything written after a space --> "copy")
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}

    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))

    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix \n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)

    # perform over-representation hypergeometic test on pathways for all patients and concatenate all results

    reactomeALL = hypergeometricTestOnMultipleIndividualsMCF(geneTable)

    write.table(reactomeALL, file=paste0("over-representationREACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)

  }
  
}
#######################################################################
#######################################################################
#######################################################################




#------------------------------- FOLDER LOOP - PATHWAY ENRICHMENT (LOGISTIC REGRESSION ON CONTINUOUS P-VALUES) -----------------------------------------------

if("-lr" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("\nPATHWAY ENRICHMENT ANALYSIS - OVER-REPRESENTATION TEST - ", folderName, "\n")
    
    # load and edit relevant samples list (remove anything written after a space --> "copy")
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
    
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix \n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    cat("loading over-representation tables \n")
    reactomeALL = read.delim(paste0("over-representationREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    keggALL = read.delim(paste0("over-representationKEGG_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    gobpALL = read.delim(paste0("over-representationGOBP_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # create PATIENTxPATHWAY matrix and fill it with adjusted p-values from the over-representation tests performed on each patient
    
    cat("creating REACTOME PATIENTxPATHWAY matrix - ", folderName, "\n")
    reactomeTable = createPatientPathwayMatrixFromHypergeometricTestOutput(reactomeALL, geneTable, sampleID, criterion)
    
    cat("creating KEGG PATIENTxPATHWAY matrix - ", folderName, "\n")
    keggTable = createPatientPathwayMatrixFromHypergeometricTestOutput(keggALL, geneTable, sampleID, criterion)
    
    cat("creating GOBP PATIENTxPATHWAY matrix - ", folderName, "\n")
    gobpTable = createPatientPathwayMatrixFromHypergeometricTestOutput(gobpALL, geneTable, sampleID, criterion)
    
    write.table(reactomeTable, file=paste0("matrixPatientPathway_LR_REACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(keggTable, file=paste0("matrixPatientPathway_LR_KEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(gobpTable, file=paste0("matrixPatientPathway_LR_GOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
    ############################
    
    # perform logistic regression with "case_ctrl" as outcome and pathways as predictors (values for each pathway are
    # adjusted p-values from hypergeometric test)
    
    cat("logistic regression on REACTOME PATIENTxPATHWAY matrix - ", folderName, "\n")
    LRmatrixREACTOME = logisticRegressionOnPatientPathwayMatrix(reactomeTable, reactomeALL, case, ctrl)
    
    cat("logistic regression on KEGG PATIENTxPATHWAY matrix - ", folderName, "\n")
    LRmatrixKEGG = logisticRegressionOnPatientPathwayMatrix(keggTable, keggALL, case, ctrl)
    
    cat("logistic regression on GOBP PATIENTxPATHWAY matrix - ", folderName, "\n")
    LRmatrixGOBP = logisticRegressionOnPatientPathwayMatrix(gobpTable, gobpALL, case, ctrl)
    
    write.table(LRmatrixREACTOME, file=paste0("logisticRegressionREACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRmatrixKEGG, file=paste0("logisticRegressionKEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRmatrixGOBP, file=paste0("logisticRegressionGOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
  
}





#######################################################################
#######################################################################
##### TEMPORARY FOR MCF!!!!!!!!!!!
#######################################################################
#------------------------------- FOLDER LOOP - PATHWAY ENRICHMENT (LOGISTIC REGRESSION ON CONTINUOUS P-VALUES) -----------------------------------------------
logisticRegressionOnPatientPathwayMatrixMCF <- function(pathwayTable, databaseALL, case, ctrl){
  
  # CHANGE CASE_CTRL (the OUTCOME) TO 0 AND 1, USE SAMPLE NAMES AS ROWNAMES, REMOVE FIRST COLUMN
  cat("pre-processing for LR\n")
  pathwayTable = pathwayTable[(pathwayTable$case_ctrl == case) | (pathwayTable$case_ctrl == ctrl),]
  pathwayTable[pathwayTable$case_ctrl == ctrl, "case_ctrl"] = 0
  pathwayTable[pathwayTable$case_ctrl == case, "case_ctrl"] = 1
  rownames(pathwayTable) = pathwayTable[,1]
  pathwayTable = pathwayTable[,-1]
  
  # perform Logistic Regression with "case_ctrl" as outcome and each pathway's adjusted p-values as predictor
  # save results in LRmatrix (matrix for results of Logistic Regression)
  
  CI=NULL #<--------------------------------------------------------------------------------------- for CI
  LRmatrix=NULL
  for(c in 2:dim(pathwayTable)[2]){
    if(c%%100 == 0){cat("fitting LR models - ", c, "/", dim(pathwayTable)[2], "\n")}
    df = pathwayTable[,c(1,c)]
    fit <- glm(case_ctrl ~ ., family = binomial(link = 'logit'), data=df)
    results = as.data.frame(coef(summary(fit)))
    t = try(confint(fit, level=0.999999))
    if(!inherits(t, "try-error")){CI = rbind(CI, as.data.frame(confint(fit, level=0.999999))[-1,])} #<-------------------------------------------- for CI
    LRmatrix = rbind(LRmatrix, results[-1,])}
  
  # adjust p-values according to BH method, sort pathwayTable according to q-values
  colnames(LRmatrix)[4] = "PVAL"
  LRmatrix$QVAL<-round(p.adjust(LRmatrix$PVAL, method="BH"), digits=20)
  # add columns for CI in LRmatrix
  LRmatrix = cbind(LRmatrix, one=rep(NA,dim(LRmatrix)[1]), two=rep(NA,dim(LRmatrix)[1]))
  colnames(LRmatrix)[(dim(LRmatrix)[2]-1):dim(LRmatrix)[2]] = c("CI_low", "CI_up")
  # copy CI values in LRmatrix, only for pathways common to both tables
  sel = which(rownames(LRmatrix) %in% rownames(CI))
  if(length(sel)>0){LRmatrix[sel, c("CI_low", "CI_up")] = CI[rownames(LRmatrix)[sel], c(1, 2)]}
  
  # add column and compute if CI comprise 1
  LRmatrix$one_in_CI_99.9999 = rep(NA, dim(LRmatrix)[1])
  for(r in 1:dim(LRmatrix)[1]){
    if(!TRUE %in% is.na(LRmatrix[r,c("CI_low", "CI_up")])){
      if(LRmatrix[r,"CI_low"]<1 & LRmatrix[r,"CI_up"]>1){LRmatrix[r,"one_in_CI_99.9999"] = "yes"
      }else{LRmatrix[r,"one_in_CI_99.9999"] = "no"}}}
  
  # sort table according to adjusted p-value
  LRmatrix<-LRmatrix[order(LRmatrix$QVAL),]
  
  # remove "`" from pathway IDs, which remains after logistic regression in REACTOME and GOBP
  if(TRUE %in% grepl("`",rownames(LRmatrix))){rownames(LRmatrix)<-unlist(lapply(strsplit(rownames(LRmatrix), "`"), "[[", 2))}
  
  # create pathwayTable with pathway IDs and pathway name "pathId2name" 
  pathId2name = unique(databaseALL[,c(1,2)])
  rownames(pathId2name) = pathId2name$ID
  
  # add columns to LRmatrix and fill with pathway IDs and assign pathway names
  # using the info in "pathId2name" pathwayTable
  cat("converting pathway IDs to pathway names\n")
  LRmatrix = cbind(pathID=rownames(LRmatrix), pathName=rep(".", dim(LRmatrix)[1]), LRmatrix, stringsAsFactors=F)
  for(r in 1:dim(LRmatrix)[1]){LRmatrix$pathName[r] = pathId2name[LRmatrix$pathID[r],2]}
  
  return(LRmatrix)
  
}

if("-lrMCF" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("\nPATHWAY ENRICHMENT ANALYSIS - OVER-REPRESENTATION TEST - ", folderName, "\n")
    
    # load and edit relevant samples list (remove anything written after a space --> "copy")
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
    
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix \n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    cat("loading over-representation tables \n")
    reactomeALL = read.delim(paste0("over-representationREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # create PATIENTxPATHWAY matrix and fill it with adjusted p-values from the over-representation tests performed on each patient
    
    cat("creating REACTOME PATIENTxPATHWAY matrix - ", folderName, "\n")
    reactomeTable = createPatientPathwayMatrixFromHypergeometricTestOutput(reactomeALL, geneTable, sampleID, criterion)
    
    write.table(reactomeTable, file=paste0("matrixPatientPathway_LR_REACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
    ############################
    
    # perform logistic regression with "case_ctrl" as outcome and pathways as predictors (values for each pathway are
    # adjusted p-values from hypergeometric test)
    
    cat("logistic regression on REACTOME PATIENTxPATHWAY matrix - ", folderName, "\n")
    LRmatrixREACTOME = logisticRegressionOnPatientPathwayMatrixMCF(reactomeTable, reactomeALL, case, ctrl)
    
    write.table(LRmatrixREACTOME, file=paste0("logisticRegressionREACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
  
}
#######################################################################
#######################################################################
#######################################################################



#------------------------------ FOLDER LOOP - PATHWAY ENRICHMENT (FISHER (1) - from OVER-REPRESENTATION - with CUT-OFF ON P-VALUES) -----------------------------------------------

if("-pf1" %in% args){
  
  # stop if no cut-off is given in command line
  if(!TRUE %in% grepl("cutoff", args)){cat("\n no cut-off specified - interrupting \n")
  }else{
    
    print(args)
    
    # extract cut-off from command line arguments
    cutoff = args[which(grepl("cutoff", args))]
    print(cutoff)
    cutoff = as.numeric(strsplit(cutoff, "=")[[1]][2])
    

    print(cutoff)
    
    for(f in 1:length(folders)){
      
      setwd(paste0(wd, "/", folders[f]))
      
      # extract relevant info from folder name
      folderName = folders[f]
      study = strsplit(folderName, "_")[[1]][1]
      platform = strsplit(folderName, "_")[[1]][2]
      criterion = strsplit(folderName, "_")[[1]][3]
      event = strsplit(folderName, "_")[[1]][4]
      ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
      case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
      
      cat("\nPATHWAY ENRICHMENT ANALYSIS - OVER-REPRESENTATION TEST - ", folderName, "\n")
      
      # load and edit relevant samples list (remove anything written after a space --> "copy")
      if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
      }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
      
      sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
      
      # load PATIENTxGENE table
      cat("loading PATIENTxGENE matrix \n")
      geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      
      cat("loading over-representation tables \n")
      reactomeALL = read.delim(paste0("over-representationREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      keggALL = read.delim(paste0("over-representationKEGG_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      gobpALL = read.delim(paste0("over-representationGOBP_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      
      
      # take adjusted p-values from the over-representation tests performed on each patient
      # create PATIENTxPATHWAY matrix and fill it with 0 if p-value >= "cutoff", or with 1 if p-value < "cutoff"
      
      cat("creating REACTOME PATIENTxPATHWAY matrix - ", folderName, "\n")
      reactomeTable = createPatientPathwayMatrixFromHypergeometricTestOutputBinary(reactomeALL, geneTable, sampleID, criterion, cutoff)
      
      cat("creating KEGG PATIENTxPATHWAY matrix - ", folderName, "\n")
      keggTable = createPatientPathwayMatrixFromHypergeometricTestOutputBinary(keggALL, geneTable, sampleID, criterion, cutoff)
      
      cat("creating GOBP PATIENTxPATHWAY matrix - ", folderName, "\n")
      gobpTable = createPatientPathwayMatrixFromHypergeometricTestOutputBinary(gobpALL, geneTable, sampleID, criterion, cutoff)
      
      write.table(reactomeTable, file=paste0("matrixPatientPathway_F1_REACTOME_", cutoff, "_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      write.table(keggTable, file=paste0("matrixPatientPathway_F1_KEGG_", cutoff, "_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      write.table(gobpTable, file=paste0("matrixPatientPathway_F1_GOBP_", cutoff, "_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      
      ############################

      # perform fisher's test on pathways between cases and ctrls
      
      cat("fisher on REACTOME PATIENTxPATHWAY matrix - ", folderName, "\n")
      RmatrixREACTOME = pathwayCaseCtrlFisher(reactomeTable, reactomeALL, case, ctrl)

      cat("fisher on KEGG PATIENTxPATHWAY matrix - ", folderName, "\n")
      RmatrixKEGG = pathwayCaseCtrlFisher(keggTable, keggALL, case, ctrl)
      
      cat("fisher on GOBP PATIENTxPATHWAY matrix - ", folderName, "\n")
      RmatrixGOBP = pathwayCaseCtrlFisher(gobpTable, gobpALL, case, ctrl)
      
      write.table(RmatrixREACTOME, file=paste0("pathwayFisher1REACTOME_", cutoff, "_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      write.table(RmatrixKEGG, file=paste0("pathwayFisher1KEGG_", cutoff, "_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      write.table(RmatrixGOBP, file=paste0("pathwayFisher1GOBP_", cutoff, "_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      
    }
    
  }
  
}


#------------------------------ FOLDER LOOP - PATHWAY ENRICHMENT (FISHER (2) - on PATHWAYS with AT LEAST ONE GENE ALTERED) -----------------------------------------------


if("-pf2" %in% args){
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    cat("   \n")
    cat("FISHER ON PATHWAYS/PATIENTS - ", folderName, "\n")
    
    # load and edit relevant samples list (remove anything written after a space --> "copy")
    cat("loading and editing sample table\n")
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
    
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    # create list of cases, ctrls and NAs according to the specific classification criterion
    selCase = separateCaseCtrl(criterion, case, ctrl, sampleID)[[1]]
    selCtrl = separateCaseCtrl(criterion, case, ctrl, sampleID)[[2]]
    selNA = separateCaseCtrl(criterion, case, ctrl, sampleID)[[3]]
    
    # load PATIENTxGENE table
    cat("loading PATIENTxGENE matrix\n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)

    ###########################

    # patient names as rownames
    rownames(geneTable) = geneTable$Sample

    # save patient names in object
    patients = geneTable$Sample

    # get all entrez IDs from gene symbols (columns)
    cat("getting entrez IDs from gene symbols\n")
    genes = colnames(geneTable)[-1]
    ids = as.vector(unlist(mget(genes,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
    ids = ids[complete.cases(ids)]

    # map all ids to reactome pathways and get unique pathway list
    cat("mapping ids to REACTOME\n")
    reactomeIds = select(reactome.db, keys=ids, column="PATHID", keytype = "ENTREZID")
    reactomeIds = unique(reactomeIds$PATHID)
    reactomeIds = reactomeIds[complete.cases(reactomeIds)]
    # get pathwayID:entrezID table for reactome
    reactome2entr = select(reactome.db, keys=reactomeIds, column="ENTREZID", keytype = "PATHID")

    # map all ids to kegg pathways and get unique pathway list
    cat("mapping ids to KEGG\n")
    keggIds = select(org.Hs.eg.db, keys=ids, column="PATH", keytype = "ENTREZID")
    keggIds = unique(keggIds$PATH)
    keggIds = keggIds[complete.cases(keggIds)]
    # get pathwayID:entrezID table for kegg
    kegg2entr = select(org.Hs.eg.db, keys=keggIds, column="ENTREZID", keytype = "PATH")

    # map all ids to gobp pathways and get unique pathway list
    cat("mapping ids to GOBP\n")
    gobpIds = select(org.Hs.eg.db, keys=ids, column="GO", keytype = "ENTREZID")
    gobpIds = gobpIds[gobpIds$ONTOLOGY == "BP",]
    gobpIds = unique(gobpIds$GO)
    gobpIds = gobpIds[complete.cases(gobpIds)]
    # get pathwayID:entrezID table for gobp
    gobp2entr = select(org.Hs.eg.db, keys=gobpIds, column="ENTREZID", keytype = "GO")

    # create PATIENTxPATHWAY matrices for REACTOME, KEGG and GOBP
    cat("creating PATIENTxPATHWAY matrix for REACTOME\n")
    pathTableREACTOME = fillPatientPathwayTable(reactomeIds, patients, reactome2entr)
    cat("creating PATIENTxPATHWAY matrix for KEGG\n")
    pathTableKEGG = fillPatientPathwayTable(keggIds, patients, kegg2entr)
    cat("creating PATIENTxPATHWAY matrix for GOBP\n")
    pathTableGOBP = fillPatientPathwayTable(gobpIds, patients, gobp2entr)

    # write.table(pathTableREACTOME, file=paste0("matrixPatientPathway_F2_REACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    # write.table(pathTableKEGG, file=paste0("matrixPatientPathway_F2_KEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    # write.table(pathTableGOBP, file=paste0("matrixPatientPathway_F2_GOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
    # if you load from file, remember to set patient names to rownames!
    # pathTableREACTOME = read.delim(paste0("matrixPatientPathwayREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    # pathTableKEGG = read.delim(paste0("matrixPatientPathwayKEGG_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    # pathTableGOBP = read.delim(paste0("matrixPatientPathwayGOBP_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # fisher's test on pathways altered in at least one gene in CASES vs CTRLS
    cat("performing fisher on REACTOME pathways\n")
    RmatrixREACTOME = pathwayCaseCtrlFisher(pathTableREACTOME, selCase, selCtrl)
    cat("performing fisher on KEGG pathways\n")
    RmatrixKEGG = pathwayCaseCtrlFisher(pathTableKEGG, selCase, selCtrl)
    cat("performing fisher on GOBP pathways\n")
    RmatrixGOBP = pathwayCaseCtrlFisher(pathTableGOBP, selCase, selCtrl)
    
    # add and edit reactome pathway names
    RmatrixREACTOME$PATHNAME = select(reactome.db, keys=RmatrixREACTOME$PATHID, column="PATHNAME", keytype ="PATHID")[,2]
    RmatrixREACTOME$PATHNAME = substr(RmatrixREACTOME$PATHNAME, 15, nchar(RmatrixREACTOME$PATHNAME))
    
    # add and edit kegg pathway names
    RmatrixKEGG[,1] = paste0("hsa", RmatrixKEGG[,1]) # adjust KEGG pathway name to be compatible to KEGGREST
    for(r in 1:dim(RmatrixKEGG)[1]){
      cat("adding KEGG pathway ID to name", r, "/", dim(RmatrixKEGG)[1], "\n")
      t = try(keggGet(RmatrixKEGG$PATHID[r])[[1]]$NAME)
      if(!inherits(t, "try-error")){RmatrixKEGG$PATHNAME[r] = keggGet(RmatrixKEGG$PATHID[r])[[1]]$NAME}} # prevent loop to stop in case no name is found
    RmatrixKEGG = RmatrixKEGG[complete.cases(RmatrixKEGG),]
    RmatrixKEGG$PATHNAME = substr(RmatrixKEGG$PATHNAME, 0, nchar(RmatrixKEGG$PATHNAME)-23)
    
    # add gobp pathway names
    RmatrixGOBP$PATHNAME = select(GO.db, keys=RmatrixGOBP$PATHID, column="TERM", keytype ="GOID")[,2]

    write.table(RmatrixREACTOME, file=paste0("pathwayPatientFisherREACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(RmatrixKEGG, file=paste0("pathwayPatientFisherKEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(RmatrixGOBP, file=paste0("pathwayPatientFisherGOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
  
}



#----------------- FOLDER LOOP - ANNOTATE GENES FROM GOterms AND COUNT HOW MANY TIMES THEY ARE ALTERED IN CASES AND CTRLS ----------------

# !!!!!!!!!!!!
# ATTENZIONE!! SE CAMBI ANALISI FISHER SUI PATHWAY, QUESTO INFLUIRÀ SU PARTE DI QUESTO SCRIPT!!!!
# !!!!!!!!!!!!

if("-a" %in% args){
  
  for(f in 1:length(folders)){
    
    setwd(paste0(wd, "/", folders[f]))
    
    # extract relevant info from folder name
    folderName = folders[f]
    study = strsplit(folderName, "_")[[1]][1]
    platform = strsplit(folderName, "_")[[1]][2]
    criterion = strsplit(folderName, "_")[[1]][3]
    event = strsplit(folderName, "_")[[1]][4]
    ctrl = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folderName, "_")[[1]][5], "vs")[[1]][2])
    
    # load PATIENTxGENE matrix
    cat("loading PATIENTxGENE matrix and sample classification table \n")
    geneTable = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # load sample classification file and edit sample names
    if(platform=="array"){sampleID <- read.delim(list.files(getwd(), 'arrayNames'), stringsAsFactors=FALSE, check.names=F)
    }else{sampleID <- read.delim("wesNames.tsv", stringsAsFactors=FALSE, check.names=F)}
    sampleID$Sample<-unlist(lapply(strsplit(sampleID$Sample, " "), "[[", 1))
    
    
    cat("\nANNOTATE GENES FOR EACH PATHWAY - LOGISTIC REGRESSION - ", folderName, "\n")

    # load and edit tables with results of logistic regression for REACTOME, KEGG, GOBP
    cat("loading logistic regression results for REACTOME \n")
    LRmatrixREACTOME = read.delim(paste0("logisticRegressionREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    cat("loading logistic regression results for KEGG \n")
    LRmatrixKEGG = read.delim(paste0("logisticRegressionKEGG_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    LRmatrixKEGG$pathID = substr(LRmatrixKEGG$pathID, 4, nchar(LRmatrixKEGG$pathID))
    cat("loading logistic regression results for GOBP \n")
    LRmatrixGOBP = read.delim(paste0("logisticRegressionGOBP_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    # create pathwayID:symbol table for REACTOME, KEGG, GOBP
    cat("getting genes contained in each pathway for REACTOME", folderName, "\n")
    reactome2symbol = select(reactome.db, keys=as.character(LRmatrixREACTOME$pathID), column="ENTREZID", keytype = "PATHID") # output is a dataframe with entrez ids in the second column
    reactome2symbol$symbol = select(org.Hs.eg.db, keys=reactome2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
    reactome2symbol = reactome2symbol[complete.cases(reactome2symbol),]
    colnames(reactome2symbol)[1] = "pathID"
    
    cat("getting genes contained in each pathway for KEGG", folderName, "\n")
    kegg2symbol = select(org.Hs.eg.db, keys=as.character(LRmatrixKEGG$pathID), column="ENTREZID", keytype = "PATH")
    kegg2symbol$symbol = select(org.Hs.eg.db, keys=kegg2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
    kegg2symbol = kegg2symbol[complete.cases(kegg2symbol),]
    colnames(kegg2symbol)[1] = "pathID"
    
    cat("getting genes contained in each pathway for GOBP", folderName, "\n")
    gobp2symbol = select(org.Hs.eg.db, keys=as.character(LRmatrixGOBP$pathID), column="ENTREZID", keytype = "GO")
    gobp2symbol = gobp2symbol[,-c(2,3)]
    gobp2symbol$symbol = select(org.Hs.eg.db, keys=gobp2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
    gobp2symbol = gobp2symbol[complete.cases(gobp2symbol),]
    colnames(gobp2symbol)[1] = "pathID"
    
    cat("expanding LRmatrix to have a row for each gene in each pathway for REACTOME", folderName, "\n")
    LRmatrixREACTOME2 = getGeneSymbolsAndCreateNewRmatrix(LRmatrixREACTOME, reactome2symbol)
    cat("expanding LRmatrix to have a row for each gene in each pathway for KEGG", folderName, "\n")
    LRmatrixKEGG2 = getGeneSymbolsAndCreateNewRmatrix(LRmatrixKEGG, kegg2symbol)
    cat("expanding LRmatrix to have a row for each gene in each pathway for GOBP", folderName, "\n")
    LRmatrixGOBP2 = getGeneSymbolsAndCreateNewRmatrix(LRmatrixGOBP, gobp2symbol)
    
    cat("calculating frequency and percentage of genes alterations in cases and ctrl for REACTOME", folderName, "\n")
    LRmatrixREACTOME2 = geneFrequencyAndPercentageInCasesAndCtrls(LRmatrixREACTOME2, geneTable, sampleID, criterion, case, ctrl)
    cat("calculating frequency and percentage of genes alterations in cases and ctrl for KEGG", folderName, "\n")
    LRmatrixKEGG2 = geneFrequencyAndPercentageInCasesAndCtrls(LRmatrixKEGG2, geneTable, sampleID, criterion, case, ctrl)
    cat("calculating frequency and percentage of genes alterations in cases and ctrl for GOBP", folderName, "\n")
    LRmatrixGOBP2 = geneFrequencyAndPercentageInCasesAndCtrls(LRmatrixGOBP2, geneTable, sampleID, criterion, case, ctrl)
    
    write.table(LRmatrixREACTOME2, file=paste0("logisticRegression_genes_REACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRmatrixKEGG2, file=paste0("logisticRegression_genes_KEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRmatrixGOBP2, file=paste0("logisticRegression_genes_GOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
    # if any file named "pathwayPatientFisher...." exist in folder
    if(length(list.files(getwd(), "pathwayPatientFisher"))>1){
      
      cat("\nANNOTATE GENES FOR EACH PATHWAY - FISHER ON PATHWAY/PATIENT - ", folderName, "\n")
      
      # load and edit tables with results of pathway/patient fisher for REACTOME, KEGG, GOBP
      cat("loading pathway/patient fisher results for REACTOME \n")
      RmatrixREACTOME = read.delim(paste0("pathwayPatientFisherREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      cat("loading pathway/patient fisher results for KEGG \n")
      RmatrixKEGG = read.delim(paste0("pathwayPatientFisherKEGG_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      RmatrixKEGG$PATHID = substr(RmatrixKEGG$PATHID, 4, nchar(RmatrixKEGG$PATHID))
      cat("loading pathway/patient fisher results for GOBP \n")
      RmatrixGOBP = read.delim(paste0("pathwayPatientFisherGOBP_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      colnames(RmatrixREACTOME)[1] = "pathID"
      colnames(RmatrixKEGG)[1] = "pathID"
      colnames(RmatrixGOBP)[1] = "pathID"
      
      # create pathwayID:symbol table for REACTOME, KEGG, GOBP
      cat("getting genes contained in each pathway for REACTOME", folderName, "\n")
      reactome2symbol = select(reactome.db, keys=as.character(RmatrixREACTOME$PATHID), column="ENTREZID", keytype = "PATHID") # output is a dataframe with entrez ids in the second column
      reactome2symbol$symbol = select(org.Hs.eg.db, keys=reactome2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
      reactome2symbol = reactome2symbol[complete.cases(reactome2symbol),]
      colnames(reactome2symbol)[1] = "pathID"
      
      cat("getting genes contained in each pathway for KEGG", folderName, "\n")
      kegg2symbol = select(org.Hs.eg.db, keys=as.character(RmatrixKEGG$PATHID), column="ENTREZID", keytype = "PATH")
      kegg2symbol$symbol = select(org.Hs.eg.db, keys=kegg2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
      kegg2symbol = kegg2symbol[complete.cases(kegg2symbol),]
      colnames(kegg2symbol)[1] = "pathID"
      
      cat("getting genes contained in each pathway for GOBP", folderName, "\n")
      gobp2symbol = select(org.Hs.eg.db, keys=as.character(RmatrixGOBP$PATHID), column="ENTREZID", keytype = "GO")
      gobp2symbol = gobp2symbol[,-c(2,3)]
      gobp2symbol$symbol = select(org.Hs.eg.db, keys=gobp2symbol$ENTREZID, column="SYMBOL", keytype = "ENTREZID")[,2]
      gobp2symbol = gobp2symbol[complete.cases(gobp2symbol),]
      colnames(gobp2symbol)[1] = "pathID"
      
      cat("expanding Rmatrix to have a row for each gene in each pathway for REACTOME", folderName, "\n")
      RmatrixREACTOME2 = getGeneSymbolsAndCreateNewRmatrix(RmatrixREACTOME, reactome2symbol)
      cat("expanding Rmatrix to have a row for each gene in each pathway for KEGG", folderName, "\n")
      RmatrixKEGG2 = getGeneSymbolsAndCreateNewRmatrix(RmatrixKEGG, kegg2symbol)
      cat("expanding Rmatrix to have a row for each gene in each pathway for GOBP", folderName, "\n")
      RmatrixGOBP2 = getGeneSymbolsAndCreateNewRmatrix(RmatrixGOBP, gobp2symbol)
      
      cat("calculating frequency and percentage of genes alterations in cases and ctrl for REACTOME", folderName, "\n")
      RmatrixREACTOME2 = geneFrequencyAndPercentageInCasesAndCtrls(RmatrixREACTOME2, geneTable, sampleID, criterion, case, ctrl)
      cat("calculating frequency and percentage of genes alterations in cases and ctrl for KEGG", folderName, "\n")
      RmatrixKEGG2 = geneFrequencyAndPercentageInCasesAndCtrls(RmatrixKEGG2, geneTable, sampleID, criterion, case, ctrl)
      cat("calculating frequency and percentage of genes alterations in cases and ctrl for GOBP", folderName, "\n")
      RmatrixGOBP2 = geneFrequencyAndPercentageInCasesAndCtrls(RmatrixGOBP2, geneTable, sampleID, criterion, case, ctrl)
      
      write.table(RmatrixREACTOME2, file=paste0("pathwayPatientFisher_genes_REACTOME_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      write.table(RmatrixKEGG2, file=paste0("pathwayPatientFisher_genes_KEGG_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      write.table(RmatrixGOBP2, file=paste0("pathwayPatientFisher_genes_GOBP_", folderName, ".tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
      
      
    }
   
       
  }
  
}



#----------------- EDIT RESULT TABLES FOR ANNA - GENE FISHER TABLES ----------------

# this option should be written like "-anna_efr criterion=x" where "x" is the criterion
# this is needed to select only the folders of interest

if("-anna_efr" %in% args){ #efr = edit fisher results
  
  # stop if no criterion is given in command line
  if(!TRUE %in% grepl("criterion", args)){cat("\n no criterion specified - interrupting \n")
  }else{
    
    # extract criterion from command line arguments
    criterion = args[which(grepl("criterion", args))]
    criterion = strsplit(criterion, "=")[[1]][2]
    
    cat("\ncriterion is", criterion, "\n")
    
    # select only folders that contain "criterion"
    folders = list.files(getwd())
    folders = folders[grep(criterion, folders)]
    
    cat("\n selected folders are \n")
    print(folders)
    
    # load and edit sample list (it's the same for all analyses, so it's ok to load it once)
    cat("\nloading sample list \n")
    first_folder = paste0(wd, "/", folders[1], "/")
    sampleID = read.delim(paste0(first_folder, list.files(first_folder, 'arrayNames')), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    sampleID$Sample = sapply(sampleID$Sample, function(x) strsplit(x, " ")[[1]][1])
    
    # get ctrl and case number info
    ctrl = as.numeric(strsplit(strsplit(folders[1], "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folders[1], "_")[[1]][5], "vs")[[1]][2])
    
    # get list of cases and ctrls
    selCtrl = sampleID[sampleID[,criterion]==ctrl, "Sample"]
    selCase = sampleID[sampleID[,criterion]==case, "Sample"]
    
    # get total number of cases and ctrls
    totCtrl = length(selCtrl)
    totCase = length(selCase)
    
    Rmatrix=NULL
    events=NULL
    
    for(f in 1:length(folders)){
      
      setwd(paste0(wd, "/", folders[f]))
      
      # extract relevant info from folder name
      folderName = folders[f]
      event = strsplit(folderName, "_")[[1]][4]
      events = append(events, event)
      
      # load gene Fisher matrix, add "event" column and concatenate to matrices for other events
      cat("concatenating gene Fisher matrix for", event, "\n")
      df = read.delim(paste0("geneFisher_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      df$event = rep(event, dim(df)[1])
      Rmatrix = rbind(Rmatrix, df)
      
      # load PATIENTxGENE matrix and assign a name according to the event
      cat("loading PATIENTxGENE matrix for", event, "\n")
      df = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      rownames(df) = df$Sample
      if(event == "gain-dup"){assign("gainDup", df)}
      else if(event == "loss-del"){assign("lossDel", df)}
      else if(event == "gain-dup-loss-del"){assign("gdld", df)}
      else{assign(event, df)}
    }
    
    # get unique gene list and its length
    genes = unique(Rmatrix$GENE)
    l = length(genes)
    
    if("gain-dup-loss-del" %in% events){events[events=="gain-dup-loss-del"] = "gdld"}
    
    # add chromosomal coordinates and cytoband info to genes
    mergedRmatrix = getCoordinatesForSymbols(as.data.frame(genes, stringsAsFactors=F), 1, entrezID="onePerSymbol")
    
    # remove "strand" info added by previous function
    mergedRmatrix = mergedRmatrix[,-which(colnames(mergedRmatrix) == "strand")]
    
    # add empty columns to summarize results for all events for all genes
    mergedRmatrix = cbind(mergedRmatrix, as.data.frame(matrix(NA, nrow=l, ncol=4*length(events))), stringsAsFactors = F)
    for(e in 1:length(events)){colnames(mergedRmatrix)[ (4*(e-1)+8) : (4*e+7) ] = c(paste0("%case_",events[e]), paste0("%ctrl_",events[e]), paste0("p_",events[e]), paste0("p-Adj_",events[e]))}
    rownames(mergedRmatrix) = genes
    
    for(r in 1:l){
      if(r%%100==0){cat("filling values for gene", r, "/", l, "\n")}
      gene = mergedRmatrix$gene[r]
      if(gene %in% colnames(gain)){mergedRmatrix[gene,8:11] = Rmatrix[Rmatrix$GENE==gene & Rmatrix$event=="gain", c("%Case", "%Ctrl", "PVAL", "QVAL")]}
      if(gene %in% colnames(gainDup)){mergedRmatrix[gene,12:15] = Rmatrix[Rmatrix$GENE==gene & Rmatrix$event=="gain-dup", c("%Case", "%Ctrl", "PVAL", "QVAL")]}
      if(gene %in% colnames(gdld)){mergedRmatrix[gene,16:19] = Rmatrix[Rmatrix$GENE==gene & Rmatrix$event=="gain-dup-loss-del", c("%Case", "%Ctrl", "PVAL", "QVAL")]}
      if(gene %in% colnames(loss)){mergedRmatrix[gene,20:23] = Rmatrix[Rmatrix$GENE==gene & Rmatrix$event=="loss", c("%Case", "%Ctrl", "PVAL", "QVAL")]}
      if(gene %in% colnames(lossDel)){mergedRmatrix[gene,24:27] = Rmatrix[Rmatrix$GENE==gene & Rmatrix$event=="loss-del", c("%Case", "%Ctrl", "PVAL", "QVAL")]}
    }
    
    write.table(mergedRmatrix, file=paste0(wd, "/geneFisherResults_anna.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
}


#----------------- EDIT RESULT TABLES FOR ANNA - PATHWAY LOG. REGR. TABLES ----------------

lookUpReactomeHigherHierarchicalLevel <- function(child, hierarchy){
  #cat("\nLOOKUP", child, "\n")
  highest=child
  
  if(length(child)>1){
    #print("if child")
    for(c in 1:length(child)){
      
      while(child[c] %in% hierarchy$child){
        
        h=NULL
        parent = hierarchy[hierarchy$child==child[c], "parent"]
        #cat("parent", parent, "\n")
        
        if(length(parent)>1){
          print("if parent")
          for(p in 1:length(parent)){
            highest = lookUpReactomeHigherHierarchicalLevel(parent[p], hierarchy)
            h = unique(append(h, highest))
            #cat("h1", h, "\n")
          }
          #cat("child1", child, "\n")
          highest = h
          #cat("highest1", highest, "\n")
        }
        
        else{
          #print("else parent")
          h = unique(append(h, parent))
          #cat("h2", h, "\n")
          highest = h
          #cat("highest2", highest, "\n")
          child = highest
          #cat("child2", child, "\n")
        }
        
      }
    }
    
  }else{
    #print("else child")
    while(child %in% hierarchy$child){
      
      h=NULL
      parent = hierarchy[hierarchy$child==child, "parent"]
      #cat("parent", parent, "\n")
      
      if(length(parent)>1){
        #print("if parent")
        for(p in 1:length(parent)){
          highest = lookUpReactomeHigherHierarchicalLevel(parent[p], hierarchy)
          h = unique(append(h, highest))
          #cat("h3", h, "\n")
        }
        #cat("child3", child, "\n")
        highest = h
        #cat("highest3", highest, "\n")
        child = NA}
      
      else{
        #print("else parent")
        h = unique(append(h, parent))
        #cat("h4", h, "\n")
        child = highest
        #cat("child4", child, "\n")
        highest = h
        #cat("highest4", highest, "\n")
      }
      
    }
  }
  
  #cat("highest5", highest, "\n\n")
  return(highest)
}

# this option should be written like "-anna_elrr criterion=x" where "x" is the criterion
# this is needed to select only the folders of interest

if("-anna_elrr" %in% args){ #efr = edit logistic regression results
  
  # stop if no criterion is given in command line
  if(!TRUE %in% grepl("criterion", args)){cat("\n no criterion specified - interrupting \n")
  }else{
    
    # extract criterion from command line arguments
    criterion = args[which(grepl("criterion", args))]
    criterion = strsplit(criterion, "=")[[1]][2]
    
    cat("\ncriterion is", criterion, "\n")
    
    # select only folders that contain "criterion"
    folders = list.files(getwd())
    folders = folders[grep(criterion, folders)]
    
    cat("\n selected folders are \n")
    print(folders)
    
    # load and edit sample list (it's the same for all analyses, so it's ok to load it once)
    cat("\nloading sample list \n\n")
    sampleID = read.delim(paste0(wd, "/", folders[1], "/arrayNames.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    sampleID$Sample = sapply(sampleID$Sample, function(x) strsplit(x, " ")[[1]][1])
    
    # get ctrl and case number info
    ctrl = as.numeric(strsplit(strsplit(folders[1], "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folders[1], "_")[[1]][5], "vs")[[1]][2])
    
    # get list of cases and ctrls
    selCtrl = sampleID[sampleID[,criterion]==ctrl, "Sample"]
    selCase = sampleID[sampleID[,criterion]==case, "Sample"]
    
    # get total number of cases and ctrls
    totCtrl = length(selCtrl)
    totCase = length(selCase)

    cat("loading REACTOME and KEGG pathway and hierarchy tables\n\n")
    # load and edit REACTOME hierarchy file
    reac_hierarchy = read.delim(paste0(pathDIR, "/reactome_hierarchy_human.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    reac_hierarchy = as.data.frame(sapply(reac_hierarchy, as.character), stringsAsFactors=F)
    # load KEGG hierarchy file
    kegg_hierarchy = read.delim(paste0(pathDIR, "/kegg_pathway_hierarchy_edited.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    # load pathID-pathName-entrezID-symbol table for REACTOME and KEGG
    reactome_path = read.delim(paste0(pathDIR, "/170428_pathID_pathName_entrezID_symbol_REACTOME.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    kegg_path = read.delim(paste0(pathDIR, "/170428_pathID_pathName_entrezID_symbol_KEGG.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    
    LRreactome=NULL
    LRkegg=NULL

    for(f in 1:length(folders)){
      
      setwd(paste0(wd, "/", folders[f]))
      
      # extract relevant info from folder name
      folderName = folders[f]
      event = strsplit(folderName, "_")[[1]][4]
      
      # load REACTOME logistic regression results, add "event" column and concatenate to matrices for other events
      cat("concatenating REACTOME Logistic Regression results for", event, "\n")
      df = read.delim(paste0("logisticRegressionREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      df$event = rep(event, dim(df)[1])
      LRreactome = rbind(LRreactome, df)

      # load KEGG logistic regression results, add "event" column and concatenate to matrices for other events
      cat("concatenating KEGG Logistic Regression results for", event, "\n")
      df = read.delim(paste0("logisticRegressionKEGG_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      df$event = rep(event, dim(df)[1])
      LRkegg = rbind(LRkegg, df)
      
      }

    # remove "R-HSA-" prefix from REACTOME pathID (new package updates display the pathID like that)
    if(TRUE %in% grepl("R-HSA-", LRreactome$pathID)){LRreactome$pathID = substr(LRreactome$pathID, 7, nchar(LRreactome$pathID))} 
    
    # remove unnecessary columns (older tables don't have CI calculations)
    LRreactome = LRreactome[,-c(4,5,6,8,9)]
    LRkegg = LRkegg[,-c(4,5,6,8,9)]
    
    # convert pathID to character (KEGG pathIDs are already character)
    LRreactome$pathID = as.character(LRreactome$pathID)
    
    LRreactome2=NULL
    LRkegg2=NULL
    
    for(f in 1:length(folders)){
      
      setwd(paste0(wd, "/", folders[f]))
      
      # extract relevant info from folder name
      folderName = folders[f]
      event = strsplit(folderName, "_")[[1]][4]
      
      # load PATIENTxGENE matrix and assign a name according to the event
      cat("\nloading PATIENTxGENE matrix for", event, "\n")
      geneFisher = read.delim(paste0("geneFisher_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)

      # get only pathways for "event"
      LRr = LRreactome[LRreactome$event == event,]
      LRk = LRkegg[LRkegg$event == event,]

      # get number of rows of LRr and LRk
      len_reac = dim(LRr)[1]
      len_kegg = dim(LRk)[1]
      
      # add empty columns to LRr
      LRr$highestLevel = rep(NA, len_reac)
      LRr$diff_altered_genes = rep(NA, len_reac)
      LRr$all_genes = rep(NA, len_reac)
      LRr$gene_proportion_1 = rep(NA, len_reac)
      LRr$gene_proportion_2 = rep(NA, len_reac)
      
      # add empty columns to LRk
      LRk$level_B = rep(NA, len_kegg)
      LRk$level_A = rep(NA, len_kegg)
      LRk$diff_altered_genes = rep(NA, len_kegg)
      LRk$all_genes = rep(NA, len_kegg)
      LRk$gene_proportion_1 = rep(NA, len_kegg)
      LRk$gene_proportion_2 = rep(NA, len_kegg)
      
      # loop over rows of LRr
      for(r in 1:len_reac ){
        
        if(r%%100==0){cat("filling info for REACTOME pathway", r, "/", len_reac, "-", event, "\n")}
        
        # get "pathway"
        pathway = LRr$pathID[r]

        # get pathID of highest levels in which "pathway" is contained
        highest = lookUpReactomeHigherHierarchicalLevel(pathway, reac_hierarchy)
        # convert highest pathID to pathname (in table "reactome_path" pathIDs and pathNames are repeated in multiple rows)
        highest = sapply(highest, function(x) unique(reactome_path[reactome_path$pathID==x, "pathName"]))
        # fill pathName of highest level
        LRr$highestLevel[r] = paste(highest, collapse = ", ")

        # get all symbols belonging to "pathway"
        symbols = reactome_path[reactome_path$pathID==pathway, "symbol"]
        # among genes in "pathway" select only those which had qval<0.0001
        alt_symbols = geneFisher[(geneFisher$GENE %in% symbols) & (geneFisher$QVAL<1e-04), "GENE"]
        
        # fill symbols of altered gene
        LRr$diff_altered_genes[r] = paste(alt_symbols, collapse=", ")
        # fill all genes
        LRr$all_genes[r] = paste(symbols, collapse = ", ")
        # fill "number of altered genes/all genes"
        LRr$gene_proportion_1[r] = paste0(length(alt_symbols), "/", length(symbols))
        # fill proportion of altered genes
        LRr$gene_proportion_2[r] = length(alt_symbols)/length(symbols)
      }
      
      # loop over rows of LRk
      for(r in 1:len_kegg){
        
        if(r%%100==0){cat("filling info for KEGG pathway", r, "/", len_kegg, "-", event, "\n")}
        
        # get "pathway"
        pathway = LRk$pathID[r]
        
        # fill name of level B
        B = kegg_hierarchy[kegg_hierarchy[, "pathID"] == pathway, "level_B"]
        if(length(B)>0){LRk$level_B[r] = B}
        # fill name of level A
        A = kegg_hierarchy[kegg_hierarchy[, "pathID"] == pathway, "level_A"]
        if(length(A)>0){LRk$level_A[r] = A}
        
        # get all symbols belonging to "pathway"
        symbols = kegg_path[kegg_path$pathID==pathway, "symbol"]
        # among genes in "pathway" select only those which had qval<0.0001
        alt_symbols = geneFisher[(geneFisher$GENE %in% symbols) & (geneFisher$QVAL<1e-04), "GENE"]
        
        # fill symbols of altered gene
        LRk$diff_altered_genes[r] = paste(alt_symbols, collapse=", ")
        # fill all genes
        LRk$all_genes[r] = paste(symbols, collapse = ", ")
        # fill "number of altered genes/all genes"
        LRk$gene_proportion_1[r] = paste0(length(alt_symbols), "/", length(symbols))
        # fill proportion of altered genes
        LRk$gene_proportion_2[r] = length(alt_symbols)/length(symbols)
        
      }
      
      # append LRr and LRk to LRreactome2 and LRkegg2
      LRreactome2 = rbind(LRreactome2, LRr)
      LRkegg2 = rbind(LRkegg2, LRk)
      
    }  
    

    # expand LRmatrix to have a row for each gene in each pathway
    LRreactome = getGeneSymbolsAndCreateNewRmatrix(LRreactome, reactome_path) # for REACTOME
    LRkegg = getGeneSymbolsAndCreateNewRmatrix(LRkegg, kegg_path)             # for KEGG


    LRreactome3=NULL
    LRkegg3=NULL

    for(f in 1:length(folders)){

      setwd(paste0(wd, "/", folders[f]))

      # extract relevant info from folder name
      folderName = folders[f]
      event = strsplit(folderName, "_")[[1]][4]

      # load PATIENTxGENE matrix and assign a name according to the event
      cat("\nloading PATIENTxGENE matrix for", event, "\n")
      matPatGene = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      # get patient names
      patients = matPatGene$Sample
      # set patient names as rownames
      rownames(matPatGene) = patients

      # for "event", compute its percentage in cases and ctrls for REACTOME
      cat("computing percentage in cases and ctrls for", event, "in REACTOME\n")
      LRr = computePercentageInCaseAndCtrlForEvent(LRreactome[LRreactome$event == event,], matPatGene, selCase, selCtrl)

      # for "event", compute its percentage in cases and ctrls for KEGG
      cat("computing percentage in cases and ctrls for", event, "in KEGG\n")
      LRk = computePercentageInCaseAndCtrlForEvent(LRkegg[LRkegg$event == event,], matPatGene, selCase, selCtrl)

      LRreactome3 = rbind(LRreactome3, LRr)
      LRkegg3 = rbind(LRkegg3, LRk)

    }

    LRreactome = LRreactome3
    LRkegg = LRkegg3
    rm(LRreactome3)
    rm(LRkegg3)

    # subset for relevant info (otherwise the file would be 500MB!!!)
    LRreactome = LRreactome[LRreactome$QVAL<0.05 & LRreactome$one_in_CI=="no", ]
    LRkegg = LRkegg[LRkegg$QVAL<0.05 & LRkegg$one_in_CI=="no", ]

    write.table(LRreactome, file=paste0(wd, "/logisticRegression_allResults_REACTOME_1_anna.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRreactome2, file=paste0(wd, "/logisticRegression_allResults_REACTOME_2_anna.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRkegg, file=paste0(wd, "/logisticRegression_allResults_KEGG_1_anna.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRkegg2, file=paste0(wd, "/logisticRegression_allResults_KEGG_2_anna.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }
 
}




#######################################################################
#######################################################################
##### TEMPORARY FOR MCF!!!!!!!!!!!
#######################################################################
#----------------- EDIT RESULT TABLES FOR GENES AND PATHWAYS - TEMPORARY FOR MCF!!!!!!!!!!!! ----------------
if("-MCF_er" %in% args){ #efr = edit logistic regression results
  
  # stop if no criterion is given in command line
  if(!TRUE %in% grepl("criterion", args)){cat("\n no criterion specified - interrupting \n")
  }else{
    
    # extract criterion from command line arguments
    criterion = args[which(grepl("criterion", args))]
    criterion = strsplit(criterion, "=")[[1]][2]
    
    cat("\ncriterion is", criterion, "\n")
    
    # select only folders that contain "criterion"
    folders = list.files(getwd())
    folders = folders[grep(criterion, folders)]
    
    cat("\n selected folders are \n")
    print(folders)
    
    # load and edit sample list (it's the same for all analyses, so it's ok to load it once)
    cat("\nloading sample list \n\n")
    sampleID = read.delim(paste0(wd, "/", folders[1], "/arrayNames.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    sampleID$Sample = sapply(sampleID$Sample, function(x) strsplit(x, " ")[[1]][1])
    
    # get ctrl and case number info
    ctrl = as.numeric(strsplit(strsplit(folders[1], "_")[[1]][5], "vs")[[1]][1])
    case = as.numeric(strsplit(strsplit(folders[1], "_")[[1]][5], "vs")[[1]][2])
    
    # get list of cases and ctrls
    selCtrl = sampleID[sampleID[,criterion]==ctrl, "Sample"]
    selCase = sampleID[sampleID[,criterion]==case, "Sample"]
    
    # get total number of cases and ctrls
    totCtrl = length(selCtrl)
    totCase = length(selCase)
    
    cat("loading REACTOME pathway and hierarchy tables\n\n")
    # load and edit REACTOME hierarchy file
    reac_hierarchy = read.delim(paste0(pathDIR, "/reactome_hierarchy_human.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    reac_hierarchy = as.data.frame(sapply(reac_hierarchy, as.character), stringsAsFactors=F)
    # load pathID-pathName-entrezID-symbol table for REACTOME and KEGG
    reactome_path = read.delim(paste0(pathDIR, "/170428_pathID_pathName_entrezID_symbol_REACTOME.tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    
    LRreactome=NULL
    
    for(f in 1:length(folders)){
      
      setwd(paste0(wd, "/", folders[f]))
      
      # extract relevant info from folder name
      folderName = folders[f]
      event = strsplit(folderName, "_")[[1]][4]
      
      # load REACTOME logistic regression results, add "event" column and concatenate to matrices for other events
      cat("concatenating REACTOME Logistic Regression results for", event, "\n")
      df = read.delim(paste0("logisticRegressionREACTOME_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      df$event = rep(event, dim(df)[1])
      LRreactome = rbind(LRreactome, df)
      
    }
    
    # remove "R-HSA-" prefix from pathID
    if(TRUE %in% grepl("R-HSA-", LRreactome$pathID)){LRreactome$pathID = substr(LRreactome$pathID, 7, nchar(LRreactome$pathID))} 
    
    # remove unnecessary columns (older tables don't have CI calculations)
    LRreactome = LRreactome[,-c(4,5,6,8,9)]
    
    # convert pathID to character (KEGG pathIDs are already character)
    LRreactome$pathID = as.character(LRreactome$pathID)
    
    LRreactome2=NULL
    
    for(f in 1:length(folders)){
      
      setwd(paste0(wd, "/", folders[f]))
      
      # extract relevant info from folder name
      folderName = folders[f]
      event = strsplit(folderName, "_")[[1]][4]
      
      # load PATIENTxGENE matrix and assign a name according to the event
      cat("\nloading PATIENTxGENE matrix for", event, "\n")
      matPatGene = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
      rownames(matPatGene) = matPatGene$Sample
      
      # get only pathways for "event"
      LRr = LRreactome[LRreactome$event == event,]
      
      # get number of rows of LRr and LRk
      len_reac = dim(LRr)[1]
      
      # add empty columns to LRr
      LRr$highestLevel = rep(NA, len_reac)
      LRr$diff_altered_genes = rep(NA, len_reac)
      LRr$all_genes = rep(NA, len_reac)
      LRr$gene_proportion_1 = rep(NA, len_reac)
      LRr$gene_proportion_2 = rep(NA, len_reac)
      
      # loop over rows of LRr
      for(r in 1:len_reac ){
        
        if(r%%100==0){cat("filling info for REACTOME pathway", r, "/", len_reac, "-", event, "\n")}
        
        # get "pathway"
        pathway = LRr$pathID[r]
        
        # get pathID of highest levels in which "pathway" is contained
        highest = lookUpReactomeHigherHierarchicalLevel(pathway, reac_hierarchy)
        # convert highest pathID to pathname (in table "reactome_path" pathIDs and pathNames are repeated in multiple rows)
        highest = sapply(highest, function(x) unique(reactome_path[reactome_path$pathID==x, "pathName"]))
        # fill pathName of highest level
        LRr$highestLevel[r] = paste(highest, collapse = ", ")
        
        # get all symbols belonging to "pathway"
        symbols = reactome_path[reactome_path$pathID==pathway, "symbol"]
        # select only genes that had a specific event among genes in "pathway"
        alt_symbols = symbols[which(symbols %in% colnames(matPatGene))]
        # fill symbols of altered gene
        LRr$diff_altered_genes[r] = paste(alt_symbols, collapse=", ")
        # fill all genes
        LRr$all_genes[r] = paste(symbols, collapse = ", ")
        # fill "number of altered genes/all genes"
        LRr$gene_proportion_1[r] = paste0(length(alt_symbols), "/", length(symbols))
        # fill proportion of altered genes
        LRr$gene_proportion_2[r] = length(alt_symbols)/length(symbols)
      }
      
      # append LRr and LRk to LRreactome2 and LRkegg2
      LRreactome2 = rbind(LRreactome2, LRr)
      
    }  
    
    
    # # expand LRmatrix to have a row for each gene in each pathway
    # LRreactome = getGeneSymbolsAndCreateNewRmatrix(LRreactome, reactome_path) # for REACTOME
    # 
    # LRreactome3=NULL
    # 
    # for(f in 1:length(folders)){
    #   
    #   setwd(paste0(wd, "/", folders[f]))
    #   
    #   # extract relevant info from folder name
    #   folderName = folders[f]
    #   event = strsplit(folderName, "_")[[1]][4]
    #   
    #   # load PATIENTxGENE matrix and assign a name according to the event
    #   cat("\nloading PATIENTxGENE matrix for", event, "\n")
    #   matPatGene = read.delim(paste0("matrixPatientGene_allOnes_", folderName,".tsv"), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
    #   # get patient names
    #   patients = matPatGene$Sample
    #   # set patient names as rownames
    #   rownames(matPatGene) = patients
    #   
    #   # for "event", compute its percentage in cases and ctrls for REACTOME
    #   cat("computing percentage in cases and ctrls for", event, "in REACTOME\n")
    #   LRr = computePercentageInCaseAndCtrlForEvent(LRreactome[LRreactome$event == event,], matPatGene, selCase, selCtrl)
    #   
    #   LRreactome3 = rbind(LRreactome3, LRr)
    # 
    # }
    # 
    # LRreactome = LRreactome3
    # rm(LRreactome3)
    # 
    # # subset for relevant info (otherwise the file would be 500MB!!!)
    # LRreactome = LRreactome[LRreactome$QVAL<0.05 & LRreactome$one_in_CI=="no", ]
    
    # write.table(LRreactome, file=paste0(wd, "/logisticRegression_allResults_REACTOME_1.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    write.table(LRreactome2, file=paste0(wd, "/logisticRegression_allResults_REACTOME_2.tsv"), sep="\t", row.names=F, col.names=T, dec=".", quote=F)
  }
  
}
#######################################################################
#######################################################################
#######################################################################

#----------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------- END OF SCRIPT ------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------

cat("\n################\nsessionInfo\n################\n")

print(sessionInfo())

cat("\n################\nend of script\n################\n")

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------- TEST -----------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------


# for protein annotation
# ids = select(org.Hs.eg.db, keys=mergedRmatrix[1:100,1], column="ENTREZID", keytype = "SYMBOL")[,2]
# ids = ids[complete.cases(ids)]
# 
# ensemble = mget(ids,envir=org.Hs.egENSEMBLPROT, ifnotfound=NA)
# 
# ncbi = mget(ids,envir=org.Hs.egREFSEQ, ifnotfound=NA)
# ncbi = sapply(prot, function(x) x = x[grepl("NP_", x) | grepl("XP_", x)])
# 
# org.Hs.egREFSEQ mergedRmatrix[1,1]


# # test version
# lookUpReactomeHigherHierarchicalLevel <- function(child, hierarchy){
#   cat("\nLOOKUP", child, "\n")
#   highest=child
#   
#   if(length(child)>1){
#     print("if child")
#     for(c in 1:length(child)){
#       
#       while(child[c] %in% hierarchy$child){
#         
#         h=NULL
#         parent = hierarchy[hierarchy$child==child[c], "parent"]
#         cat("parent", parent, "\n")
#         
#         if(length(parent)>1){
#           print("if parent")
#           for(p in 1:length(parent)){
#             highest = lookUpReactomeHigherHierarchicalLevel(parent[p], hierarchy)
#             h = unique(append(h, highest))
#             cat("h1", h, "\n")
#           }
#           cat("child1", child, "\n")
#           highest = h
#           cat("highest1", highest, "\n")
#         }
#         
#         else{
#           print("else parent")
#           h = unique(append(h, parent))
#           cat("h2", h, "\n")
#           highest = h
#           cat("highest2", highest, "\n")
#           child = highest
#           cat("child2", child, "\n")
#         }
#         
#       }
#     }
#     
#   }else{
#     print("else child")
#     while(child %in% hierarchy$child){
#       
#       h=NULL
#       parent = hierarchy[hierarchy$child==child, "parent"]
#       cat("parent", parent, "\n")
#       
#       if(length(parent)>1){
#         print("if parent")
#         for(p in 1:length(parent)){
#           highest = lookUpReactomeHigherHierarchicalLevel(parent[p], hierarchy)
#           h = unique(append(h, highest))
#           cat("h3", h, "\n")
#         }
#         cat("child3", child, "\n")
#         highest = h
#         cat("highest3", highest, "\n")
#         child = NA}
#       
#       else{
#         print("else parent")
#         h = unique(append(h, parent))
#         cat("h4", h, "\n")
#         child = highest
#         cat("child4", child, "\n")
#         highest = h
#         cat("highest4", highest, "\n")
#       }
#       
#     }
#   }
#   
#   cat("highest5", highest, "\n\n")
#   return(highest)
# }
# 
# top = lookUpReactomeHigherHierarchicalLevel(c("177929","422475"), reac_hierarchy)
# top = sapply(top, function(x) all_path[all_path$PATHID==x, "PATHNAME"])


# setwd("/home/PERSONALE/eugenio.fonzi2/53M/amltp53_array_53M_gain-dup-loss-del_1vs2")
#setwd("/Users/Eugenio/Desktop/Lavoro/Seragnoli/anna/tp53/test/170412_53M_all")
# 
# del = read.delim(paste0("arrayTable_del.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
# dup = read.delim(paste0("arrayTable_dup.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)
# gainDup = read.delim("amltp53_array_53M_gain-dup_1vs2/matrixPatientGene_allOnes_amltp53_array_53M_gain-dup_1vs2.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F)
# gdld = read.delim("amltp53_array_53M_gain-dup-loss-del_1vs2/matrixPatientGene_allOnes_amltp53_array_53M_gain-dup-loss-del_1vs2.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F)
# loss2 = read.delim("loss.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F)
# matrix = read.delim(paste0("matrixPatientGene_allOnes_amltp53_array_53M_gain-dup-loss-del_1vs2.tsv"), header=T, sep="\t", stringsAsFactors = F, check.names = F)



#FALSE %in% (c("a","b","c") == c("a","b","c"))
