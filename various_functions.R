hypergeometricTestOnMultipleIndividuals_01 <- function(geneTable){

  # loop over patients and perform over-representation test for REACTOME, KEGG, GOBP pathways (based on hypergeometric distribution)
  # requires a PATIENTxGENE matrix of zeroes and ones "geneTable"

  reactomeALL=NULL
  keggALL=NULL
  gobpALL=NULL
  
  # loop over table rows (patients)
  for(r in 1:dim(geneTable)[1]){
    
    # get name of the patient
    patient = geneTable[r,1]
    
    # select altered genes for the patient
    genes = colnames(geneTable)[which(geneTable[r,] == 1)]
    
    cat("     \n")
    cat("PATIENT", r, "of", dim(geneTable)[1], " - ", length(genes), "events", "\n")
    
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
      
      # perform over-representation test for kegg, add sample name, add case/ctrl classification number and
      # concatenate to "keggALL"
      if(length(ids)>0){
        kegg = try(enrichKEGG(ids, organism="hsa", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1))
        if(!inherits(kegg, "NULL")){
          kegg = kegg@result
          kegg$Sample = patient
          keggALL = rbind(keggALL, kegg)
          cat("kegg --> DONE \n")
        }else{cat("SKIP kegg \n")}}
      
      # perform over-representation test for gobp, add sample name, add case/ctrl classification number and
      # concatenate to "gobpALL"
      # try() and inherits() allow to decide whether to skip the rest of the loop in case enrichGO() finds an error and returns NULL
      if(length(ids)>0){
        gobp = try(enrichGO(gene=ids, OrgDb=org.Hs.eg.db, keytype = "ENTREZID", ont = "BP", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1))
        if(!inherits(gobp, "NULL")){
          gobp = gobp@result
          gobp$Sample = patient
          gobpALL = rbind(gobpALL, gobp)
          cat("gobp --> DONE \n")
        }else{cat("SKIP gobp\n")}}
      
    }
  }
  
  return(list(reactomeALL, keggALL, gobpALL))
  
}
