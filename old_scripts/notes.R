bamList = system("ls wes_anna/batches/double_germ/Analysis/Data-cleanup/BQSR/ | grep .bam", intern = T)

for(b in 1:length(bamList)){  system(paste0("> ~/wes_anna/trial10/", bamList[b]))}

varList = unique(vars$SAMPLE)
for(v in 1:length(varList)){system(paste0("> ~/wes_anna/trial10/", varList[v]))}

args = c("~/wes_anna/trial7", 
         "~/wes_anna/sample_pairing.tsv", 
         "~/reference", 
         "~/tools/jdk1.6.0_35/jre/bin/java", 
         "~/tools/mutect/muTect-1.1.4.jar", 
         "~/tools/VarScan.v2.3.9.jar",
         "~/tools/samtools_0.1.2/samtools")
p=2

kit="truseq"

######## FOR GENE-PATHWAY ENRICHMENT

wd = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/aneuploidy/"
wd = '/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/test/170605_aneupl_compl_karyo'
setwd(wd)

folders = list.files(getwd())

sampleID = read.delim("170606_arrayNames.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F)

sampleID2 = read.delim("170606_arrayNames.tsv", header=T, sep="\t", stringsAsFactors = F, check.names = F)

sampleID3 = read.delim(list.files(getwd(), 'arrayNames'), header=T, sep="\t", stringsAsFactors = F, check.names = F)


first_folder = paste0(wd, "/", folders[1], "/")
sampleID4 = read.delim(paste0(first_folder, list.files(first_folder, 'arrayNames')), header=T, sep="\t", dec=".", stringsAsFactors = F, check.names=F)
