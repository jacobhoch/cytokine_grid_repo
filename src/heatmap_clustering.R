# kmeans based on naive cyto DEGs

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')

args          <- parser$parse_args()
filestem      <- args$s

set.seed(123)
dds_list <- readRDS(file = glue("./data/DE_results/{filestem}_list.Rds"))
vsd <- assays(dds_list[[1]])[['vsd']]
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

cytokines <- getCytokines()
cyto_DEGs <- list()
for (i in cytokines) {
  # all DEGs
  #DE_info <- read.csv(file = glue("~/data/DE_results/{filestem}/none_{i}_vs_none_none_DE.csv"))
  #cyto_DEGs[[i]] <- DE_info$X
  
  # only up DEGs
  cyto_DEGs[[i]] <- readtxt(glue("./data/DE_results/{filestem}/none_{i}_vs_none_none_up.txt"))
}

combined_DEGs <- Reduce(union, cyto_DEGs)

# cluster by sample type
unique_DEGs <- list()
for (i in cytokines) {
  unique_DEGs[[i]] <- cyto_DEGs[[i]][!(cyto_DEGs[[i]] %in%
                          unlist(cyto_DEGs[which(names(cyto_DEGs) != i)], use.names=F))]
}
plot_vsd <- vsd[rownames(vsd) %in% unlist(unique_DEGs, use.names=F),]
z <- t(scale(t(plot_vsd)))
gene_label <- enframe(unique_DEGs) %>%
  unnest(cols = value) %>%
  as.data.frame()
rownames(gene_label) <- gene_label$value
gene_label <- gene_label[order(match(rownames(gene_label),rownames(z))),]
gene_label <- gene_label$name

pdf(file=glue("./fig/sfigs/cluster_by_cytokine_DEGs.pdf"),
    width = 9, height = 8)
ha <- HeatmapAnnotation(cytokine1 = sub(".*:(.*?)_.*", "\\1", colnames(plot_vsd)),
                        cytokine2 = sub(".*_","", colnames(plot_vsd)),
                        col = list(cytokine1 = getCytoColors(),
                                   cytokine2 = getCytoColors()))
ht <- Heatmap(z, top_annotation = ha, show_row_names = F, row_split = gene_label,
              show_column_names = F, column_km = 5, row_title_rot = 0,
              column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
ht <- draw(ht)
dev.off()
