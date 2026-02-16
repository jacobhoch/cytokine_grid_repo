# now looking at concentration data

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-f', '--filestem', type='character', nargs=1, help='stem where DEG results live')

args          <- parser$parse_args()
filestem      <- args$f

safe_extract <- function(lst, n) {
  if (length(lst) >= n) {
    return(lst[[n]])
  } else {
    return(NA)
  }
}

concColor <- function() {col <- c("H"="#CD87F8","M"="#E3BCFB","L"="#F5E7FE")}
timeColor <- function() {col <- c("24"="#CD87F8","6"="#E3BCFB","1"="#F5E7FE")}

# IL4 DEG heatmap --------------------------------------------------------------
set.seed(123)
dds_list <- readRDS(file=glue("./data/DE_results/{filestem}_list.Rds"))
dds <- dds_list[[1]]

vsd <- vst(dds, blind=FALSE, fitType = 'local')
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

IL4_DEG_list <- read.csv(file=glue("./data/DE_results/{filestem}/none.IL4_vs_none.none_DE.csv"))
IFNb_DEGs <- readtxt(file=glue("./data/DE_results/{filestem}/none.IFNb_H_24_vs_none.none_up.txt"))
#IFNb_DEGs <- readtxt(file=glue("./data/DE_results/{filestem}/IFNb_H_24.IL4_vs_none.none_up.txt"))
IL4_DEGs <- setdiff(IL4_DEG_list[IL4_DEG_list$log2FoldChange > 1.4, ]$X, IFNb_DEGs)

condition_list <- c("none.none","none.IL4","IFNb_H_[0-9]+.IL4",
                    "IFNb_M_[0-9]+.IL4","IFNb_L_[0-9]+.IL4")
sample_list <- colnames(vsd)[grepl(paste0(condition_list, collapse="|"),colnames(vsd))]

plot_vsd <- vsd[IL4_DEGs, sample_list]
z <- t(scale(t(plot_vsd)))

cytokine1_info <- strsplit(sub(".*\\.(.*)\\..*", "\\1", colnames(z)), split = "_")
cytokine1 <- sapply(cytokine1_info, "[[", 1)
conc1 <- unlist(lapply(cytokine1_info, safe_extract, n=2))
time1 <- unlist(lapply(cytokine1_info, safe_extract, n=3))

stim_label <- data.frame(name=substr(colnames(z),3,nchar(colnames(z))),
                         row.names=colnames(z))
stim_label$name <- factor(stim_label$name,
                          levels = c("none.none","none.IL4",
                                     "IFNb_L_1.IL4","IFNb_L_6.IL4","IFNb_L_24.IL4",
                                     "IFNb_M_1.IL4","IFNb_M_6.IL4","IFNb_M_24.IL4",
                                     "IFNb_H_1.IL4","IFNb_H_6.IL4","IFNb_H_24.IL4"))

pdf(file=glue("./fig/time_conc/IL4_DEGs_IFNb_conds.pdf"),
    width = 6, height = 9)
ha <- HeatmapAnnotation(cytokine1 = sapply(cytokine1_info, "[[", 1),
                        concentration = conc1, time = time1,
                        cytokine2 = sub("^[^.]*\\.[^.]*\\.", "\\1", colnames(z)),
                        col = list(cytokine1 = getCytoColors(),
                                   concentration = concColor(),
                                   time = timeColor(),
                                   cytokine2 = getCytoColors()),
                        simple_anno_size = unit(3, "mm"))
ht <- Heatmap(z, top_annotation = ha, column_split = stim_label, column_title = NULL,
              show_row_names = T, row_km = 2, row_title = NULL,
              show_column_names = F, column_gap = unit(0.4, "mm"), row_gap = unit(1, "mm"),
              row_names_gp = gpar(fontsize = 7.5), cluster_column_slices = FALSE,
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in'),
                                          legend_direction = "horizontal"))
ht <- draw(ht, heatmap_legend_side = "bottom")
dev.off()

### FC heatmap -----------------------------------------------------------------
# condition_list <- c("IFNb_H_24.IL4_vs_none.IFNb_H_24","IFNb_M_24.IL4_vs_none.IFNb_M_24","IFNb_L_24.IL4_vs_none.IFNb_L_24",
#                     "IFNb_H_6.IL4_vs_IFNb_H_6.none","IFNb_M_6.IL4_vs_IFNb_M_6.none","IFNb_L_6.IL4_vs_IFNb_L_6.none",
#                     "IFNb_H_1.IL4_vs_IFNb_H_1.none","IFNb_M_1.IL4_vs_IFNb_M_1.none","IFNb_L_1.IL4_vs_IFNb_L_1.none",
#                     "none.IL4_vs_none.none")
# 
# genes <- read.csv(file = glue("~/data/DE_results/{filestem}/{condition_list[1]}_full.csv"))$X
# FCmatrix <- data.frame(matrix(ncol = 0, nrow = length(genes), 
#                               dimnames = list(genes, NULL)))
# for (i in condition_list) {
#   
#   full_csv <- read.csv(file = glue("~/data/DE_results/{filestem}/{i}_full.csv"))
#   FCmatrix[[i]] <- full_csv$log2FoldChange
# }
# 
# IL4_DEG_list <- read.csv(file=glue("~/data/DE_results/{filestem}/none.IL4_vs_none.none_DE.csv"))
# IL4_DEGs <- IL4_DEG_list[IL4_DEG_list$log2FoldChange > 1.5, ]$X
# 
# FCmatrix <- FCmatrix[rownames(FCmatrix) %in% IL4_DEGs, ]
# 
# stim_label <- data.frame(substr(colnames(z),3,nchar(colnames(z))),
#                          row.names=colnames(z))
# 
# ht <- Heatmap(FCmatrix, show_row_names = T, #column_split = stim_label
#               show_column_names = F, row_km = 3, column_title_gp = gpar(fontsize=0),
#               column_gap = unit(1, "mm"), row_gap = unit(1, "mm"), row_dend_side = "left",
#               heatmap_legend_param = list(title="z score",
#                                           title_gp=gpar(fontsize=12),
#                                           labels_gp=gpar(fontsize=10),
#                                           legend_height=unit(1.5,'in'),
#                                           legend_width=unit(0.5,'in')))
# draw(ht)
