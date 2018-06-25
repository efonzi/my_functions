##############################################################################################################################
#################################################### FUNCTIONS ###############################################################
##############################################################################################################################

# compute number of months passed between pairs of dates
# it returns a 1-column matrix with all the computed months
# 'diag' contains the earlier dates, while 'lfu' contains the later dates
# both must be 3-column dataframes, with days, months and years labelled as 'dd', 'mm', 'yy' (in any order)
# ----------------
# | mm | dd | yy |
# ----------------
# |  2 | 23 | 10 |
# | 12 |  4 | 98 |
# ----------------
computeMonthsFromDate <- function(diag, lfu){
  
  # check if the two tables have the same dimensions
  if(FALSE %in% (dim(diag) == dim(lfu))){stop("diagnosis and last_follow_up tables have different dimensions!!!")}
  
  # get number of samples
  s = dim(diag)[1]
  
  # create sx1 matrix of zeroes
  months = matrix(0, s, ncol=1)
  
  # calculate no. survival months for each sample and store it in a corrispondent row in "surv_months"
  for(i in 1:s){# loop over samples (rows)
    
    # get years difference
    ydiff = lfu[i,'yy'] - diag[i,'yy']
    
    if(ydiff<0){ydiff = (100+lfu[i,'yy']-diag[i,'yy'])} # add 100 in case diagnosis date was before 2000 and FU date was after 2000
    
    # get months difference
    if(ydiff==0){mdiff = lfu[i,'mm'] - diag[i,'mm']
    }else{
      if(diag[i,'mm']>lfu[i,'mm']){mdiff = (ydiff-1)*12 + (12-diag[i,'mm']) + lfu[i,'mm']}
      else{mdiff = ydiff*12 + (lfu[i,'mm'] - diag[i,'mm'])}
    }
    
    # add value to matrix
    months[i,1] = mdiff
  }
  
  return(months)
}




# annotate gene coordinates for ANY table containing gene symbols
# "gene_col" is the number of the column containing the gene symbols
# possible arguments for "entrezID" are: 
#  > "all" in order to get all entrezID for each symbol
#  > "onePerSymbol" to get only the first entrezID for each symbol
getCoordinatesForSymbols <- function(genes, gene_col, entrezID){
  
  library(org.Hs.eg.db)
  
  # get gene list
  symbol = genes[,gene_col]
  
  # get symbol2entrezID table
  genes2entrez = select(org.Hs.eg.db, keys=as.character(symbol), columns="ENTREZID", keytype = "SYMBOL")
  colnames(genes2entrez)[2] = "entrez"
  
  if(entrezID=="all"){
    # as it has multiple entrezID per symbol, get those symbols
    sel = genes2entrez[which(duplicated(genes2entrez$SYMBOL)), "SYMBOL"]
    dup = genes[genes[,gene_col] %in% sel,]
    # then add rows for duplicated symbols to "genes"
    genes = rbind(genes, dup)
    
    # sort both tables by symbol and merge
    genes = genes[order(genes[,gene_col]),]
    genes2entrez = genes2entrez[order(genes2entrez$SYMBOL),]
    genes = cbind(genes, entrez=genes2entrez$entrez, stringsAsFactors=F)
  }
  
  if(entrezID=="onePerSymbol"){
    genes2entrez = genes2entrez[-which(duplicated(genes2entrez$SYMBOL)),]
    genes = cbind(genes, entrez=genes2entrez$entrez, stringsAsFactors=F)
  }
  
  # get all entrezID
  entrez = genes$entrez
  entrez = entrez[complete.cases(entrez)]
  
  # add columns
  genes$chr = rep(NA, dim(genes)[1])
  genes$start = rep(NA, dim(genes)[1])
  genes$end = rep(NA, dim(genes)[1])
  genes$strand = rep("+", dim(genes)[1])
  genes$chrArm = rep(NA, dim(genes)[1])
  genes$cytoband = rep(NA, dim(genes)[1])
  
  # get start positions of entrez ids
  startPos = mget(entrez,envir=org.Hs.egCHRLOC, ifnotfound=NA)
  start = sapply(startPos, function(x) x[[1]])
  start = start[!is.na(start)]
  
  # get chromosome of entrez ids
  chromosome = sapply(startPos, function(x) names(x)[1])
  chromosome = unlist(chromosome)
  
  # get end positions of entrez ids
  endPos = mget(entrez,envir=org.Hs.egCHRLOCEND, ifnotfound=NA)
  end = sapply(endPos, function(x) x[[1]])
  end = end[!is.na(end)]
  
  # get cytobands of entrez ids
  cytoband = mget(entrez,envir=org.Hs.egMAP, ifnotfound=NA)
  cyto = sapply(cytoband, function(x) x[[1]])
  cyto = cyto[!is.na(cyto)]
  
  # loop over table rows and fill start positions, end positions and strand sign of genes
  for(r in 1:dim(genes)[1]){
    if(r%%1000==0){cat("filling info for gene", r, "/", dim(genes)[1], "\n")}
    entrez = genes$entrez[r]
    if(!is.na(entrez)){
      if(entrez %in% names(chromosome)){genes$chr[r] = chromosome[entrez]}
      if(entrez %in% names(start)){genes$start[r] = as.character(start[entrez])}
      if(entrez %in% names(end)){genes$end[r] = as.character(end[entrez])}
      if(entrez %in% names(cyto)){
        genes$chrArm[r] = cyto[entrez]
        genes$cytoband[r] = cyto[entrez]}
      if(grepl("-", genes$start[r])){
        genes$start[r] = as.numeric(strsplit(genes$start[r], "-")[[1]][2]) # positions on negative strands come as negative integers
        genes$end[r] = as.numeric(strsplit(genes$end[r], "-")[[1]][2]) # positions on negative strands come as negative integers
        genes$strand[r] = "-"}
    }
  }
  
  # convert "weird" chromosome names like "1_KI270762v1_alt" to just chromosome number
  genes$chr = sapply(genes$chr, function(x) if(grepl("_",x)){strsplit(x, "_")[[1]][1]}else{x=x})
  
  # extract chromosome arm info
  genes[grep("q", genes$chrArm), "chrArm"] = as.vector(sapply(genes[grep("q", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "q")[[1]][1], "q")))
  genes[grep("p", genes$chrArm), "chrArm"] = as.vector(sapply(genes[grep("p", genes$chrArm), "chrArm"], function(x) paste0(strsplit(x, "p")[[1]][1], "p")))
  
  genes$chr = as.character(genes$chr)
  genes$start = as.numeric(genes$start)
  genes$end = as.numeric(genes$end)
  
  return(genes)
  
}





# create "segments" table with the coordinates of all altered segments
# "table" must be an edited output from NEXUS, "order" must be a table specifing the radial position of each patient and
# "palette" must be a table assigning color to specific events
createSegmentsForHighlight <- function(table, order, palette){
  
  l = dim(table)[1]
  
  highlight = data.frame(chr=rep(".",l), start=rep(0, l), end=rep(0, l), options=rep(".", l), stringsAsFactors = F)
  
  gainDup = NULL
  lossDel = NULL
  
  for(r in 1:l){
    
    # get event type
    event = table$Event[r]
    
    # get coordinates from table in the format "chr#:x,xxx,xxx-x,xxx,xxx"
    coordinates = table$Chromosome.Region[r]
    
    # extract chr number and save as "hs#"
    chr = strsplit(coordinates, ":")[[1]][1]
    chr = paste0("hs", strsplit(chr, "chr")[[1]][2])
    
    # extract "start" and "end" positions and convert to format "xxxxxxx"
    positions = strsplit(coordinates, ":")[[1]][2]
    
    start = strsplit(positions, "-")[[1]][1]
    start = as.numeric(gsub(",", "", start))
    
    end = strsplit(positions, "-")[[1]][2]
    end = as.numeric(gsub(",", "", end))
    
    # determine color of event
    color = palette[palette$event==event, "color"]
    
    # get radius coordinates (r0 and r1) for patient
    patient = table$Sample[r]
    r0 = order[order$Sample==patient, "r0"]
    r1 = order[order$Sample==patient, "r1"]
    
    # fill info for "highlight" dataset --> chromosome, start, end and options(color, r0, r1)
    highlight$chr[r] = chr
    highlight$start[r] = start
    highlight$end[r] = end
    highlight$options[r] = paste0("fill_color=", color, ",r0=", r0, "r,r1=", r1, "r")
    
    # fill info for GRanges of a specific event
    if(event == "gain" | event == "dup"){gainDup = rbind(gainDup, data.frame(chr=chr, start=start, end=end, stringsAsFactors = F))}
    if(event == "loss" | event == "del"){lossDel = rbind(lossDel, data.frame(chr=chr, start=start, end=end, stringsAsFactors = F))}
    
  }
  
  gainDup = makeGRangesFromDataFrame(gainDup)
  lossDel = makeGRangesFromDataFrame(lossDel)
  
  return(list(highlight, gainDup, lossDel))
  
}





# create co-occurrence table for LINKS --> a table with 1:1 gene combinations with counts of patients that harbor them
# PATIENTxGENE matrix for gainDup and for lossDel are necessary (any event could be suitable)
# column with gene names in "genes" must be labelled "symbol"
computeCombinationsForLinks <- function(gainDup, lossDel, genes){
  
  # get only desired genes (columns)
  links_gd = gainDup[,which(colnames(gainDup) %in% genes$symbol)]
  links_ld = lossDel[,which(colnames(lossDel) %in% genes$symbol)]
  
  # merge gain-dup table to loss-del table
  links_gdld = cbind(links_gd, links_ld)
  
  # substitute events with corresponding gene name, in order to compute combinations
  for(c in 1:dim(links_gdld)[2]){
    gene = colnames(links_gdld)[c]
    links_gdld[links_gdld[,c]>0,c] = gene}
  
  # for each patient calculate all 1:1 combinations of any type of event and create a total list
  combinations=NULL
  for(r in 1:dim(links_gdld)[1]){
    events = as.character(links_gdld[r,which(links_gdld[r,]!=0)])
    if(length(events)>1){
      df = as.data.frame(t(combn(events, 2)), stringsAsFactors=F) # combn() calculates all possible combinations of 2 events; t() transposes the resulting matrix
      combinations = rbind(combinations, df)}}
  
  colnames(combinations) = c("a","b")
  
  # count how many times each combination is repeated
  combinations = ddply(combinations,.(a,b),nrow)
  
  return(combinations) 
}






# create dataframes for links in A and links in E, with following structure:
# chr1 - start1 - end1 - chr2 - start2 - end2 - color=....
createLinksTable <- function(combinations, genes, paletteName, paletteTot, paletteType, max){
  
  # create paletteString
  paletteString = paste0("color=", paletteName, "-", paletteTot, "-", paletteType, "-")
  
  # get total number of combinations
  l = dim(combinations)[1]
  
  # create empty dataframe with necessary fields
  links = data.frame(chr1=rep(".",l), start1=rep(0,l), end1=rep(0,l), chr2=rep(".",l), start2=rep(0,l), end2=rep(0,l), options=rep(0,l), stringsAsFactors = F)
  
  # loop over rows and fill each cell
  for(r in 1:l){
    gene1 = combinations$a[r]
    gene2 = combinations$b[r]
    count = combinations$V1[r]
    links$chr1[r] = genes[genes$symbol==gene1, "chr"]
    links$start1[r] = genes[genes$symbol==gene1, "start"]
    links$end1[r] = genes[genes$symbol==gene1, "end2"]
    links$chr2[r] = genes[genes$symbol==gene2, "chr"]
    links$start2[r] = genes[genes$symbol==gene2, "start"]
    links$end2[r] = genes[genes$symbol==gene2, "end2"]
    
    # # assign color according to level of link "connection"
    # if(count < max){links$options[r] = paste0(paletteString, ((paletteTot/max)*count)%/%1+1)
    # }else{links$color[r] = paste0(paletteString, paletteTot)}}
    
    # to group the less "connected" links and assign them the same color (which is not the clearest one and it will be visible)
    if(count < max){num = ((paletteTot/max)*count)%/%1+1
    if(num > 2){links$options[r] = paste0(paletteString, num, ",z=", num)
    }else{links$options[r] = paste0(paletteString, "2,z=2")}
    }else{links$options[r] = paste0(paletteString, paletteTot, ",z=", paletteTot)}}
  
  # # to convert clear color coefficients to high color and viceversa (suitable only for scales from 1 to 11)
  # if(count < max){num = ((paletteTot/max)*count)%/%1+1
  #   links$options[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))
  # }else{num = paletteTot
  #   links$options[r] = paste0(paletteString, num+10*(1-0.2*(num-1)))}}
  
  return(links)
  
}
