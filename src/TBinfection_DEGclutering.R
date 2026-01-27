# goal is: cateogrization of genes
# methods: kmeans or just venn/upset diagrams
# trying kmeans first


source("~/src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')
parser$add_argument('-u', '--uninfected_data', type='character', nargs=1,
                    help='stem where uninfected DEG results live')

args          <- parser$parse_args()
filestem      <- args$s
uninf_stem    <- args$u

dds_list <- readRDS(glue("~/data/DE_results/{filestem}_list.Rds"))
vsd <- assays(dds_list[[1]])[['vsd']]
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

sample_list <- colnames(vsd)[grepl(glue("_{cyto}"),colnames(vsd)) |
                               grepl("none_none", colnames(vsd))]

conditions <- c("none_IL4","WT_none","EccCa1.Tn_none")
DEGs_by_condition <- list()
for (i in conditions) {
  DE_info <- read.csv(file = glue("~/data/DE_results/{filestem}/{i}_vs_none_none_DE.csv"))
  
  DEGs_by_condition[[i]] <- DE_info$X
  #cyto_DEGs[[i]] <- readtxt(glue("~/data/DE_results/{filestem}/none_{i}_vs_none_none_up.txt"))
}

combined_DEGs <- Reduce(union, DEGs_by_condition)

sample_list <- colnames(vsd)[!grepl("IFNb",colnames(vsd))]

plot_vsd <- vsd[rownames(vsd) %in% combined_DEGs,
                colnames(vsd) %in% sample_list]
z <- t(scale(t(plot_vsd)))
stim_label <- data.frame(substr(colnames(z),3,nchar(colnames(z))),
                         row.names=colnames(z))
stim_label[,1] <- sub("_"," to ",stim_label[,1])

# kmeans by gene
ha <- HeatmapAnnotation(day1 = sub(".*_(.+)_.*", "\\1", colnames(plot_vsd)),
                        day2 = sub(".*_","", colnames(plot_vsd)),
                        col = list(day1 = getInfColors(),
                                   day2 = getInfColors()))
ht <- Heatmap(z, column_split = stim_label, top_annotation = ha, show_row_names = F,
              show_column_names = F, row_km = 5, column_title_gp = gpar(fontsize=0),
              column_gap = unit(0.25, "mm"), row_gap = unit(1, "mm"),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
draw(ht)

#kmeans by sample
IL4_DEGs <- DEGs_by_condition$none_IL4
WT_DEGs <- DEGs_by_condition$WT_none

unique_DEGs <- list("IL4"=setdiff(IL4_DEGs, WT_DEGs))#,
                    #"Mtb"=setdiff(WT_DEGs, IL4_DEGs))

plot_vsd <- vsd[rownames(vsd) %in% unlist(unique_DEGs, use.names=F),
                colnames(vsd) %in% sample_list]

z <- t(scale(t(plot_vsd)))
gene_label <- enframe(unique_DEGs) %>%
  unnest(cols = value) %>%
  as.data.frame()
rownames(gene_label) <- gene_label$value
gene_label <- gene_label[order(match(rownames(gene_label),rownames(z))),]
gene_label <- gene_label$name

set.seed(123)
ht <- Heatmap(z, top_annotation = ha, show_row_names = F, row_km = 1,
              show_column_names = F, column_km = 1, row_title_rot = 0,
              column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
ht <- draw(ht)

#IL4 DEGs only, kmeans by gene
plot_vsd <- vsd[rownames(vsd) %in% IL4_DEGs,
                colnames(vsd) %in% sample_list]
z <- t(scale(t(plot_vsd)))
stim_label <- data.frame(substr(colnames(z),3,nchar(colnames(z))),
                         row.names=colnames(z))
stim_label[,1] <- sub("_"," to ",stim_label[,1])

# kmeans by gene
ha <- HeatmapAnnotation(day1 = sub(".*_(.+)_.*", "\\1", colnames(plot_vsd)),
                        day2 = sub(".*_","", colnames(plot_vsd)),
                        col = list(day1 = getInfColors(),
                                   day2 = getInfColors()))
ht <- Heatmap(z, column_split = stim_label, top_annotation = ha, show_row_names = F,
              show_column_names = F, row_km = 3, column_title_gp = gpar(fontsize=0),
              column_gap = unit(0.25, "mm"), row_gap = unit(1, "mm"),
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
draw(ht)


### VENN DIAGRAM TIME
MtbIL4_DEGs <- read.csv(file = glue("~/data/DE_results/{filestem}/WT_IL4_vs_WT_none_DE.csv"))
MtbIL4_DEGs <- MtbIL4_DEGs$X

IL4_DEGs <- read.csv(file = glue("~/data/DE_results/{uninf_stem}/none_IL4_vs_none_none_DE.csv"))
IL4_DEGs <- IL4_DEGs$X

IFNbIL4_DEGs <- read.csv(file = glue("~/data/DE_results/{uninf_stem}/IFNb_IL4_vs_IFNb_none_DE.csv"))
IFNbIL4_DEGs <- IFNbIL4_DEGs$X

colors <- c("IFNb\nto IL4"="#CD87F8","Mtb\nto IL4"="#B8D455","IL4"="#929000")
venn_genes <- list("IFNb\nto IL4"=IFNbIL4_DEGs, "Mtb\nto IL4"=MtbIL4_DEGs, "IL4"=IL4_DEGs)

venn_genes <- list("IFNb to IL4" = setdiff(IFNbIL4_DEGs, setdiff(IL4_DEGs, WT_DEGs)),
                   "Mtb to IL4" = setdiff(MtbIL4_DEGs, setdiff(IL4_DEGs, WT_DEGs)))

plot(euler(venn_genes, shape='circle'),
     fill=list(fill=colors[names(venn_genes)], alpha=0.3),
     edges=list(col=colors[names(venn_genes)], lwd=4),
     quantities=list(cex=1.8, col="black", font=2),
     labels=list(cex=2, col=colors[names(venn_genes)]))

### IFNB + WT VENN DIAGRAM TIME ------------------------------------------------
Mtb_DEGs <- read.csv(file = glue("~/data/DE_results/{filestem}/WT_none_vs_none_none_DE.csv"))
Mtb_DEGs <- Mtb_DEGs$X

IFNb_DEGs <- read.csv(file = glue("~/data/DE_results/{uninf_stem}/IFNb_none_vs_none_none_DE.csv"))
IFNb_DEGs <- IFNb_DEGs$X

EccCa1_DEGs <- read.csv(file = glue("~/data/DE_results/{filestem}/EccCa1.Tn_none_vs_none_none_DE.csv"))
EccCa1_DEGs <- EccCa1_DEGs$X

colors <- c("IFNb"="#CD87F8","Mtb"="#B8D455","EccCa1:Tn"="#31AFF5")
venn_genes <- list("IFNb"=IFNb_DEGs, "Mtb"=Mtb_DEGs, "EccCa1:Tn"=EccCa1_DEGs)

plot(euler(venn_genes, shape='circle'),
     fill=list(fill=colors[names(venn_genes)], alpha=0.3),
     edges=list(col=colors[names(venn_genes)], lwd=4),
     quantities=list(cex=1.8, col="black", font=2),
     labels=list(cex=2, col=colors[names(venn_genes)]))
