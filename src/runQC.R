plotQC <- function(d, basename, ext='pdf') {
  fig1 <- plot_total_counts(d)
  fig2 <- plot_library_complexity(d)
  fig3 <- plot_gene_detection(d)
  
  plots <- list(fig1, fig2, fig3)
  
  p <- cowplot::plot_grid(plotlist = plots, ncol = length(plots))
  cowplot::save_plot(file=glue('{basename}.{ext}'), plot=p, base_asp = length(plots))
}

source('./src/utilities.R')
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-i', '--input_file', type='character', nargs=1, help='DESeq2 dataset object')
parser$add_argument('-r', '--replicates-drop', type='character', nargs='+', default='none',
                    help='names of replicates to drop')
parser$add_argument('-d', '--design', type='character', nargs='+', default='none',
                    help='DESeq design, default is design of raw dds object')

args        <- parser$parse_args()
file        <- args$i
drop        <- args$r
alt_design  <- args$d
stem <- sub("\\..*", "",sub(".*/","",file))

dds <- readRDS(file)

# collapsing replicates
dds <- collapseReplicates(dds, groupby=dds$Condition, run=dds$Unique)

fig_stem <- './fig/QC/'
fig_name <- basename(file_path_sans_ext(file))

plotQC(dds, glue('{fig_stem}{fig_name}_pre'))

sample_threshold = ncol(dds) / 2
count_threshold = 3
# take counts data -> cpm
cpm <- assay(dds) / (colSums(assay(dds)) / 1e6)
# require gene counts to be > 3
cpm_data_threshold = rowSums(cpm > count_threshold, na.rm=T) >= sample_threshold

min_size <- getMinSizes()[stem]
size_threshold = colSums(counts(dds)) > min_size

dds_filt <- dds[cpm_data_threshold, size_threshold]
print(paste0("removed: ", colnames(dds)[!size_threshold]))

if (length(drop) > 1 & !('none' %in% drop)) {
  dds_filt <- dds_filt[, !colnames(dds_filt) %in% drop]
}

assays(dds_filt)[['vsd']] <- vst(dds_filt, fitType = "local")

plotQC(dds_filt, glue('{fig_stem}{fig_name}_post'))

# if alternate design, save different dds object
if (length(alt_design) > 1 & !('none' %in% alt_design)) {
  
  design(dds_filt) <-  as.formula(paste0("~",paste(alt_design, collapse="+")))
  design_ext <- gsub("[^[:alnum:]]", "", design(dds_filt))[2]
  
} else {
  design_ext <- ""
}

saveRDS(dds_filt, file=glue('./data/clean_dds/{fig_name}{design_ext}.Rds'))