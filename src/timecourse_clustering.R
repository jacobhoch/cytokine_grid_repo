# kmeans timepoint data

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-f', '--filestem', type='character', nargs=1, help='stem where DEG results live')

args          <- parser$parse_args()
filestem      <- args$f

dds_list <- readRDS(file=glue("./data/DE_results/{filestem}_list.Rds"))
vsd <- assays(dds_list[[1]])[['vsd']]
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)
vsd <- assay(vsd)

# IFNb DEG heatmap -------------------------------------------------------------
condition_list <- c("none.IFNb_H","none.none")
donor_list <- c("G\\.","H\\.","I\\.")
sample_list <- colnames(vsd)[grepl(paste0(condition_list, collapse="|"),
                                   colnames(vsd))]
sample_list <- sample_list[grepl(paste0(donor_list, collapse="|"), sample_list)]

IFNb24_DEGs <- readtxt(file = glue("./data/DE_results/{filestem}/none.IFNb_H_24_vs_none.none_up.txt"))
IFNb6_DEGs <- readtxt(file = glue("./data/DE_results/{filestem}/none.IFNb_H_6_vs_none.none_up.txt"))
IFNb1_DEGs <- readtxt(file = glue("./data/DE_results/{filestem}/none.IFNb_H_1_vs_none.none_up.txt"))

IFNb_DEGs <- Reduce(union, list(IFNb24_DEGs, IFNb6_DEGs, IFNb1_DEGs))

plot_vsd <- vsd[rownames(vsd) %in% IFNb_DEGs,
                colnames(vsd) %in% sample_list]
z <- t(scale(t(plot_vsd)))

pdf(file=glue("./fig/sfigs/IFNb_high_DEGs.pdf"),
    height = 10, width = 6)
stim_label <- data.frame(substr(colnames(z),3,nchar(colnames(z))),
                         row.names=colnames(z))
ht <- Heatmap(z, column_split = stim_label, show_row_names = F,
              show_column_names = F, row_km = 4, column_title_gp = gpar(fontsize=0),
              column_gap = unit(0.25, "mm"), row_gap = unit(1, "mm"), row_dend_side = "left",
              heatmap_legend_param = list(title="z score",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(1.5,'in'),
                                          legend_width=unit(0.5,'in')))
draw(ht)
dev.off()

### KMEANS CLUSTER TIMECOURSE GRAPH --------------------------------------------
KMeans <- kmeans(z, 4, iter.max = 10, nstart = 50)
Groups <- as.factor(KMeans$cluster)
samplesize <- dim(z)[1]
id <- factor(1:samplesize)

PlotCluster <- data.frame(c(z[,grepl("none.none",colnames(z))]),
                          c(z[,grepl("1",colnames(z))]),
                          c(z[,grepl("6",colnames(z))]), 
                          c(z[,grepl("24",colnames(z))]),
                          as.factor(c(id,id,id)),
                          as.factor(Groups))

colnames(PlotCluster)<-c('0','1','6','24','ID','Cluster')
PlotAllClusters<-melt(PlotCluster,by='ID')

my_colors <- c("#F4A261","#E76F51","#264653","#2A9D8F","#E9C46A") 
#my_colors <- c("#950d97","#e802db","#fd12ff","#ffbbff")

lineplot <- ggplot(PlotAllClusters,aes(variable,value,group=ID,col=Cluster)) +
  geom_point(position = "jitter", size=1,alpha=.5) +  # add jittered points with low opacity
  #geom_line(position = "jitter", alpha=.1) +
  stat_smooth(aes(group=Cluster), colour="black", 
              method="loess", se=F, size=1, formula = 'y~x') +
  ylab('z score') +
  xlab('IFNb duration (hours)') +
  facet_grid(. ~ factor(Cluster, levels=c("1","2","3","4")), 
             labeller = as_labeller(c("1" = "Cluster 1","2" = "Cluster 2",
                                      "3" = "Cluster 3","4" = "Cluster 4"))) +
  scale_color_manual(values = my_colors, 
                     breaks =c("1","2","3","4")) +
  guides(col = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        strip.text = element_blank())

cluster_summary <- data.frame(size = PlotCluster$Cluster)
cluster_summary$size <- factor(cluster_summary$size, levels=c(1, 2, 3, 4))

bar_plot <- ggplot(cluster_summary, aes(y=size, fill=size)) +
  geom_bar() +
  scale_fill_manual(values = my_colors) +
  xlab("number of genes") +
  ylab("cluster") +
  theme_classic() +
  coord_flip() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")

pdf(file = glue("./fig/time_conc/IFNb_high_DEGs_clustering.pdf"),
    width = 7, height = 5)
lineplot / bar_plot
dev.off()
