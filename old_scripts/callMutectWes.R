
# Run like: Rscript script_name.R argument1 argument2
# argument1 = folder containing BAM files written with "/" at the end; there must be a BAM for NORMAL sample and a BAM for TUMOR sample
# argument2 = path to table with normal-tumor pairings (must have columns called "normal", "tumor", "kit")
# argument3 = 
# argument4 = 
# argument5 = 

# set depending directories
refDIR = "~/reference/"
javaDIR = "~/tools/jdk1.6.0_35/jre/bin/java"
mutectDIR = "~/tools/mutect/muTect-1.1.4.jar"

# get and print arguments
args <- commandArgs(trailingOnly = TRUE)
cat("\n arguments are\n")
print(args)

# set working directory
setwd(args[1])

# get sample pairing table
pairings = read.delim(args[2], header=T, sep="\t", quote="", stringsAsFactors=F)

# get and print names of all BAM files
bamList = list.files(pattern = ".bam$")
cat("\nBAM files are\n")
print(bamList)

# get and print names of all BAI files
baiList = list.files(pattern = ".bai$")
cat("\nBAI files are\n")
print(baiList)

# BAM and BAI file names are in the form: "9_tumor_Sample_78540_sequence_1.fastq_Sample_78540_sequence_2.fastq_aligned.sam.bam"
# in which the sample ID is just "78540", as in the sample pairing file (argument2)
# there is therefore need to change file names

# # loop over pairs of sample names
# for(p in 1:dim(pairings)[1]){
#   
#   # if both normal and tumor sample names are contained (as substring) among BAM file names
#   if(TRUE %in% (grepl(pairings$normal[p], bamList)) & TRUE %in% (grepl(pairings$tumor[p], bamList))){
#     
#     # get BAM file name for normal and tumor
#     normal = bamList[which(grepl(pairings$normal[p], bamList))]
#     tumor = bamList[which(grepl(pairings$tumor[p], bamList))]
# 
#     # change BAM file names to the name written in the sample pairing table
#     system(paste0("mv ", args[1], normal, " ", args[1], pairings$normal[p], ".bam"))
#     system(paste0("mv ", args[1], tumor, " ", args[1], pairings$tumor[p], ".bam"))
#     
#     # get BAI file name for normal and tumor
#     normal = baiList[which(grepl(pairings$normal[p], baiList))]
#     tumor = baiList[which(grepl(pairings$tumor[p], baiList))]
#     
#     # change BAI file names to the name written in the sample pairing table
#     system(paste0("mv ", args[1], normal, " ", args[1], pairings$normal[p], ".bai"))
#     system(paste0("mv ", args[1], tumor, " ", args[1], pairings$tumor[p], ".bai"))
#     
#   }
#   # get names of BAM files after changing names
#   bamList = list.files(pattern = ".bam$")
#   # get names of BAI files after changing names
#   baiList = list.files(pattern = ".bai$")
# }
# 
# # print changed BAM names
# cat("\nBAM file names were changed to\n")
# print(bamList)
# 
# # print changed BAI names
# cat("\nBAI file names were changed to\n")
# print(baiList)

# loop over pair of sample names
for(p in 1:dim(pairings)[1]){

  # if both normal and tumor sample names are contained (as substring) among BAM file names
  if(TRUE %in% (grepl(pairings$normal[p], bamList)) & TRUE %in% (grepl(pairings$tumor[p], bamList))){

    # get normal and tumor sample name (no extension)
    normal = pairings$normal[p]
    tumor = pairings$tumor[p]
    # get amplification kit name
    if(pairings$kit[p]=="truseq"){kit="truseq-exome-targeted-regions-manifest-v1-2_edited.bed"}else{kit="nexterarapidcapture_expandedexome_targetedregions_fixed.bed"}

    cat("\n Normal:", normal, " - Tumor:", tumor, " - Kit:", pairings$kit[p], "\n", sep = "")

    print("Running MuTect on normal vs tumor")
    # call MuTect
    system(paste0(javaDIR, " -Xmx2g -jar ", mutectDIR, " --analysis_type MuTect --reference_sequence ", refDIR, "human_g1k_v37.fasta --cosmic ", refDIR, "b37_cosmic_v54_120711.vcf --dbsnp ", refDIR, "dbsnp_138.b37.vcf --intervals ", refDIR, kit, " --input_file:normal ", args[1], normal, ".bam --input_file:tumor ", args[1], tumor, ".bam --out ", args[1], normal, "_", tumor, ".mutect.call_stats.txt --coverage_file ", args[1], normal, "_", tumor, ".mutect.coverage.wig.txt 2>&1 > ", args[1], normal, "_", tumor, ".mutect.err_log.txt"))

  }
}

# get err_logs names
errList = list.files(pattern = ".mutect.err_log.txt$")

# concatenate err_logs
system(paste0("cat ", paste(errList, collapse = " "), " > mutect_err_log.txt"))

# remove single err_logs
system(paste0("rm ", paste(errList, collapse = " ")))


