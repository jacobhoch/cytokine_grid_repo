# IL4 DEGs wtih IFNb stim

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-f', '--filestem', type='character', nargs='*', default='1,2', help='')

args        <- parser$parse_args()
filestem    <- args$f

cyto_colors <- getCytoColors()

set.seed(123)
dds_list <- readRDS(file=glue("./data/DE_results/{filestem}_list.Rds"))
dds <- dds_list[[1]]

# saving intermediates from DESeq PCA
vsd <- vst(dds, blind=FALSE, fitType = 'local')
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

IL4_DEG_list <- read.csv(file=glue("./data/DE_results/{filestem}/none_IL4_vs_none_none_DE.csv"))
IL4_DEGs <- IL4_DEG_list[IL4_DEG_list$log2FoldChange > 1.4, ]$X

condition_list <- c("none_none","none_IL4","IFNb_none","IFNb_IL4","IL4_IL4","IL4_IFNb")
sample_list <- colnames(vsd)[grepl(paste0(condition_list, collapse="|"),colnames(vsd))]

plot_vsd <- vsd[IL4_DEGs, sample_list]
z <- t(scale(t(plot_vsd)))

pdf(file=glue("./fig/cytokine_grid/heatmap_IL4_DEGs_wIFNb.pdf"),
    width = 6, height = 9)
ha <- HeatmapAnnotation(cytokine1 = sub(".*:(.*?)_.*", "\\1", colnames(plot_vsd)),
                        cytokine2 = sub(".*_","", colnames(plot_vsd)),
                        col = list(cytokine1 = getCytoColors(),
                                   cytokine2 = getCytoColors()))
ht <- Heatmap(z, top_annotation = ha, show_row_names = T, row_km = 3, row_title = NULL,
              show_column_names = F, column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
              row_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
ht <- draw(ht)
dev.off()
