args = commandArgs()
WORKDIR=(paste0(args[6],'/'))
SAMDIR=(paste0(args[7],'/'))
setwd(WORKDIR)

filenames = list.files(pattern="bam")
for(i in 1:length(filenames)){
  system("echo '--------------------'")
  system(paste0("echo 'flagstat ", filenames[i], "'"))
  system("echo '--------------------'")
  system(paste0(SAMDIR, 'samtools flagstat ', filenames[i]))
  

}
