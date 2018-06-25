source("https://bioconductor.org/biocLite.R")

biocLite()

biocLite("BiocUpgrade") # to upgrade bioconductor version

library(BiocInstaller)

pack = "org.Hs.eg.db"
pack = "reactome.db"
pack = "clusterProfiler"
pack = "ReactomePA"
pack = "GO.db"
pack = "KEGGREST"
pack = "SomaticSignatures"
pack = "BSgenome.Hsapiens.UCSC.hg19"

biocLite(pack)

library(ReactomePA)

biocValid("clusterProfiler")
biocValid()

# functions in "GenomicRanges" can't
# from https://support.bioconductor.org/p/93347/#93373
pkgs <- c(
  "S4Vectors", "IRanges", "GenomicRanges", "DelayedArray",
  "XVector", "GenomicAlignments", "ShortRead",
  "VariantAnnotation", "AnnotationHub", "GGtools",
  "ggbio", "ChemmineR", "InteractionSet", "flowCore",
  "GenomicTuples", "CNEr", "MultiAssayExperiment",
  "genomeIntervals", "TFBSTools", "IWTomics", "spliceSites",
  "podkat", "kebabs", "matter", "dada2",
  "ClusterSignificance", "gespeR", "HiTC", "tigre", "soGGi"
)
update <- intersect(rownames(installed.packages()), pkgs)
biocLite(update, type="source", ask=FALSE)
