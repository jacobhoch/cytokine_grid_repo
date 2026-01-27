# now looking at concentration data

source("~/src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')

args          <- parser$parse_args()
filestem      <- args$s
uninf_stem    <- args$u

dds_list <- readRDS(file=glue("~/data/DE_results/{filestem}_list.Rds"))
vsd <- assays(dds_list[[1]])[['vsd']]
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

# IFNb DEG heatmap -------------------------------------------------------------
sample_list <- colnames(vsd)[grepl("IFNb_H_24.IL4",colnames(vsd)) |
                               grepl("IFNb_M_24.IL4",colnames(vsd)) |
                               grepl("IFNb_L_24.IL4",colnames(vsd)) |
                               grepl("none.IL4",colnames(vsd)) |
                               grepl("none.none", colnames(vsd)) &
                               !grepl("none_", colnames(vsd))]

IFNbH_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/IFNb_H_24.IL4_vs_none.IFNb_H_24_up.txt"))
IFNbM_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/IFNb_M_24.IL4_vs_none.IFNb_M_24_up.txt"))
IFNbL_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/IFNb_L_24.IL4_vs_none.IFNb_L_24_up.txt"))

IFNbH_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/IFNb_H_24.IL4_vs_none.IFNb_H_24_up.txt"))
IFNbM_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/IFNb_M_24.IL4_vs_none.IFNb_M_24_up.txt"))
IFNbL_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/IFNb_L_24.IL4_vs_none.IFNb_L_24_up.txt"))

IFNb_DEGs <- Reduce(union, list(IFNbH_DEGs, IFNbM_DEGs, IFNbL_DEGs))
IL4_DEGs <- readtxt(file = glue("~/data/DE_results/{filestem}/none.IL4_vs_none.none_up.txt"))

plot_vsd <- vsd[rownames(vsd) %in% IL4_DEGs,
                colnames(vsd) %in% sample_list]
z <- t(scale(t(plot_vsd)))
z <- z[, !colnames(z) %in% c("A.none_IL4","B.none_IL4","C.none_IL4")]
stim_label <- data.frame(substr(colnames(z),3,nchar(colnames(z))),
                         row.names=colnames(z))

ht <- Heatmap(z, column_split = stim_label, show_row_names = F,
              show_column_names = T, row_km = 2, column_title_gp = gpar(fontsize=0),
              column_gap = unit(1, "mm"), row_gap = unit(1, "mm"), row_dend_side = "right",
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
draw(ht)
