prepEnv <- function() {
  libraries <- c("GenomeInfoDb", "EnsDb.Hsapiens.v75", "ggplot2", "patchwork",
                 "biovizBase", "DESeq2", "edgeR", "limma", "tidyverse", "dplyr",
                 "apeglm","org.Hs.eg.db", "DOSE", "AnnotationHub", "ensembldb",
                 "biomaRt", "zeallot", "ggforce","RColorBrewer","ComplexHeatmap",
                 "glue","tools","argparse", "RNAseqQC", "EnhancedVolcano",
                 "gridExtra","UpSetR")
  
  invisible(lapply(libraries, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
  
  #complex heatmap options
  ht_opt("merge_legends" = TRUE) 
}

getCytoColors <- function() {
  cyto_colors <- c("IFNb"="#CD87F8","IFNg"="#FFD966","IL10"="#005493","IL1b"="#FF7E79",
                   "IL4"="#929000","none"="#CCCCCC","TGFb"="#941100")
}

getInfColors <- function() {
  colors <- c("none"="#CCCCCC","IFNb"="#CD87F8","IL4"="#929000","WT"="#b8d455",
              "EccCa1.Tn"="#31AFF5")
}

getCytokines <- function() {
  cytokines <- c("IFNb","IFNg","IL10","IL1b","IL4","TGFb")
}

getMinSizes <- function() {
  min_sizes <- list("241120_primeseq" = 250000, "240722_primeseq" = 150000, "250201_TBinfection" = 500000)
}

getMainRef <- function() {
  main_refs <- list("241120_primeseq" = c("Stim","none_none"),
                    "240722_primeseq" = c("Stim","none.none"),
                    "250201_TBinfection" = c("Stim","none_none"))
}

readtxt <- function(file) {
  con <- file(file, open="r")
  txt <- readLines(con)
  close(con)
  return(txt)
}

getNormalizedCounts <- function(dds) {
  counts <- counts(dds, normalized=TRUE)
}
