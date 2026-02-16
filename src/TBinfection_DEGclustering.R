source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-f', '--filestem', type='character', nargs=1, help='stem where DEG results live')

args          <- parser$parse_args()
filestem      <- args$f

dds_list <- readRDS(glue("./data/DE_results/{filestem}_list.Rds"))
vsd <- assays(dds_list[[1]])[['vsd']]
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

#kmeans by sample
IL4_DEG_list <- read.csv(file=glue("./data/DE_results/{filestem}/none_IL4_vs_none_none_DE.csv"))
IL4_DEGs <- IL4_DEG_list[IL4_DEG_list$log2FoldChange > 1.6, ]$X

sample_list <- colnames(vsd)[!grepl("IFNb",colnames(vsd))]
plot_vsd <- vsd[rownames(vsd) %in% IL4_DEGs,
                colnames(vsd) %in% sample_list]

z <- t(scale(t(plot_vsd)))

# kmeans by gene
pdf(file = glue("./fig/TBinfection/IL4_DEGs_heatmap.pdf"),
    width = 6, height = 10)
ha <- HeatmapAnnotation(infection = sub(".*_(.+)_.*", "\\1", colnames(plot_vsd)),
                        cytokine = sub(".*_","", colnames(plot_vsd)),
                        col = list(infection = getInfColors(),
                                   cytokine = getInfColors()))
ht <- Heatmap(z, top_annotation = ha, show_row_names = T, row_names_gp = gpar(fontsize=8),
              show_column_names = F, row_km = 3, row_title_gp = gpar(fontsize=0),
              column_gap = unit(0.25, "mm"), row_gap = unit(1, "mm"),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
draw(ht)
dev.off()
