# generate buckets of DEGs
# then maybe make a heatmap out of them? who's to say

getDEcategory <- function(naive_genes, repol_DEGs) {
  row <- c()
  row[1] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE", 1])
  row[2] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE likely", 1])
  row[3] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE not excluded", 1])
  row[4] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE unlikely", 1])
  return(row)
}

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')
parser$add_argument('-a', '--a_cyto', type='character', nargs=1, help='first cytokine')
parser$add_argument('-b', '--b_cyto', type='character', nargs=1, help='second cytokine')

args        <- parser$parse_args()
filestem    <- args$s
cyto_A      <- args$a
cyto_B      <- args$b

repol_DEGs <- read.csv(glue("~/data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_full.csv"))
naive_DEGs <- read.csv(glue("~/data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_full.csv"))
dds_list <- readRDS(glue("~/data/DE_results/{filestem}_list.Rds"))

vsd <- assays(dds_list[[1]])[['vsd']]
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

colnames(repol_DEGs)[10] <- "category_repol"
colnames(naive_DEGs)[10] <- "category_naive"
DEG_data <- cbind(naive_DEGs, repol_DEGs["category_repol"])

context_dep_genes <- DEG_data[DEG_data$category_naive == "DE" & DEG_data$category_repol != "DE",]
gene_list <- context_dep_genes$X

sample_list <- colnames(vsd)[grepl(glue("{cyto_A}_{cyto_B}"),colnames(vsd)) |
                               grepl(glue("none_{cyto_B}"),colnames(vsd)) |
                               grepl(glue("{cyto_A}_none"),colnames(vsd)) |
                               grepl("none_none", colnames(vsd))]

plot_vsd <- vsd[rownames(vsd) %in% gene_list,
                grepl(paste(sample_list, collapse="|"), colnames(vsd))]

z <- t(scale(t(plot_vsd)))
stim_label <- data.frame(substr(colnames(z),3,nchar(colnames(z))),
                         row.names=colnames(z))
stim_label[,1] <- sub("_"," to ",stim_label[,1])
stim_label[,1] <- factor(stim_label[,1], levels=c("none to none",
                                                  glue("none to {cyto_B}"),
                                                  glue("{cyto_A} to {cyto_B}"),
                                                  glue("{cyto_A} to none")))

png(filename = glue("~/fig/{filestem}/heatmap_{cyto_A}_to_{cyto_B}_CDGs.png"))
ht <- Heatmap(z, column_split = stim_label,
              show_column_names = F,
              column_title_side = "bottom", column_title_rot = 20,
              row_names_gp = gpar(fontsize=min(500/dim(z)[1],9)),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
draw(ht)
dev.off()

# DRAW TABLES ------------------------------------------------------------------
categories <- c("DE","DE likely","DE not excluded","DE unlikely")
grid <- data.frame(DE = rep(NA,4),
                   DE_likely = rep(NA,4),
                   DE_notexclude = rep(NA,4),
                   DE_unlikely = rep(NA,4),
                   row.names = categories)
colnames(grid) <- gsub("_"," ",colnames(grid))

for (i in categories){
  # pull gene information to determine which DEGs to look at
  naive_genes <- naive_DEGs[naive_DEGs$category == i,1]
  
  grid[i,] <- getDEcategory(naive_genes, repol_DEGs)

}
png(filename = glue("./fig/{filestem}/DEG_{cyto_A}_to_{cyto_B}_truthtable.png"),
    width = 2400, height = 700, res = 200)
tb <- tableGrob(grid, theme=ttheme_default(base_size=25))
ggplot(mapping = aes(grid)) +
  annotation_custom(tb) +
  ylab(glue("{cyto_B} vs none")) +
  xlab(glue("{cyto_A} to {cyto_B} vs {cyto_A}")) +
  theme_classic() +
  scale_x_discrete(position = "top") +
  theme(axis.title.x.top = element_text(size = 30),
        axis.title.y = element_text(size = 30))
dev.off()

