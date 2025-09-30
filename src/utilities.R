prepEnv <- function() {
  libraries <- c("GenomeInfoDb", "EnsDb.Hsapiens.v75", "ggplot2", "patchwork",
                 "biovizBase", "DESeq2", "edgeR", "limma", "tidyverse", "dplyr",
                 "apeglm","org.Hs.eg.db", "DOSE", "AnnotationHub", "ensembldb",
                 "biomaRt", "zeallot", "ggforce","RColorBrewer","ComplexHeatmap",
                 "glue","tools","argparse", "RNAseqQC")
  
  suppressPackageStartupMessages(lapply(libraries, require, character.only=TRUE))
}

getCytoColors <- function() {
  cyto_colors <- c("IFNb"="#CD87F8","IFNg"="#FFD966","IL10"="#005493","IL1b"="#FF7E79",
                   "IL4"="#929000","none"="#CCCCCC","TGFb"="#941100")
}