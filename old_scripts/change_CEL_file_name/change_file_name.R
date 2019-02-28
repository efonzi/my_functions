# args[1] must be the working directory
# args[2] must be the name of the file with the name conversion table
# args[3] must be the desired extension (like ".txt")

args <- commandArgs(trailingOnly = T)
wd <- args[1]
setwd(wd)

# load conversion table
conversion = read.delim(args[2], header=T, sep="\t", stringsAsFactors=F, check.names=F)
# remove anything written after a space
conversion = as.data.frame(apply(conversion, 2, function(x) sapply(x, function(y) strsplit(y, " ")[[1]][1])), stringsAsFactors=F)



# get names of all files with the desired extension (args[3]) and cut it out
files = sapply(list.files(getwd(), pattern=args[3]), function(x) strsplit(x, args[3])[[1]][1])
#print(head(files))

# make an output directory
system(paste0("mkdir ", args[1], "/output"))

# loop over rows in conversion files, change old file name to new file name and move file to "output" directory
if(!FALSE %in% (conversion[,"old_name"] %in% files)){# check if all needed files are present
  for(n in 1:dim(conversion)[1]){
    old = conversion[n,"old_name"]


    if(grepl("[(]", old){gsub( ##### PROBLEM WITH "(" IN FILE NAMES!!!! ---> IT SHOULD BE "\(" IN BASH!!!!!


    new = conversion[n,"new_name"]
    system(paste0("mmv ", args[1], "/", old, args[3], " ", args[1], "/", new, args[3]))
    system(paste0("mv ", args[1], "/", new, args[3], " ", args[1], "/output"))
  }
}else{# print error message if any file is missing
  print("the following files are missing:")
  print(conversion[which(!conversion[,1] %in% files),1])
}
  
