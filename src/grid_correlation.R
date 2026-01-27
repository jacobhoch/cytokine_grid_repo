# grid 1.0 clustering heatmap

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-d', '--dds_list', type='character', nargs=1, help='')

args    <- parser$parse_args()
file    <- args$d

cytokines <- getCytokines()
cyto_colors <- getCytoColors()

set.seed(123)
dds_list <- readRDS(file)
dds <- dds_list[[1]]
vsd <- vst(dds, blind=FALSE, fitType = 'local')
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)

sampleDistMatrix <- cor(assay(vsd), method = 'pearson')

colors <- colorRampPalette(rev(brewer.pal(9, "RdGy")))(255)
cyto_colors2 = c("IFNb"="#CD87F8","IFNg"="gray90","IL10"="gray82","IL1b"="gray75",
                 "IL4"="gray67","TGFb"="gray60","none"="gray50")
cyto_colors <- getCytoColors()

ha_df <- data.frame(cbind(cytokine1=as.character(vsd$day_1),
                          cytokine2=as.character(vsd$day_2)))
ha <- HeatmapAnnotation(df=ha_df,
                        col=list(cytokine1 = cyto_colors,
                                 cytokine2 = cyto_colors),
                        which = "row",
                        simple_anno_size = unit(0.5, "cm"))

kclus <- kmeans(sampleDistMatrix, 2)
split <- factor(kclus$cluster,
                levels=c("2","1"))

pdf(file=glue("./fig/cytokine_grid/correlation_heatmap.pdf"))
ht <- Heatmap(sampleDistMatrix, show_column_names=F, row_gap=unit(1.5,'mm'),
              column_gap=unit(0.1,'mm'), show_row_names=F, show_row_dend = F,
              right_annotation = ha, col=colors, use_raster = F, row_split = split,
              heatmap_legend_param = list(title = "Pearson correlation",
                                          labels_gp = gpar(fontsize=12),
                                          title_gp = gpar(fontsize=12, fontface='bold')),
              show_heatmap_legend=T, width=unit(11.5,'cm'), height=unit(11.5,'cm'),
              clustering_method_columns = "average",
              clustering_method_rows = "average")
draw(ht)
dev.off()

