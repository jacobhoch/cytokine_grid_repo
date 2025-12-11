# upset plot of all new DEG categories given previous category
source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')
parser$add_argument('-c', '--cyto', type='character', nargs=1, help='main cytokine')

args        <- parser$parse_args()
filestem    <- args$s
cyto        <- args$c

naive_gene_info <- read.csv(glue("~/data/DE_results/{filestem}/none_{cyto}_vs_none_none_full.csv"))

naive_DEGs <- naive_gene_info[naive_gene_info$category == "DE",]$X

cytokines <- getCytokines()

CDG_by_cytokine <- list()
CIG_by_cytokine <- list()
for (i in cytokines[cytokines != cyto]) {
  cyto_gene_info <- read.csv(glue("~/data/DE_results/{filestem}/{i}_{cyto}_vs_{i}_none_full.csv"))
  
  CIGs <- intersect(naive_DEGs, cyto_gene_info[cyto_gene_info$category == "DE",]$X)
  CDGs <- intersect(naive_DEGs, cyto_gene_info[cyto_gene_info$category != "DE",]$X)
  
  CIG_by_cytokine[[i]] <- CIGs
  CDG_by_cytokine[[i]] <- CDGs
}

#png(filename = glue("./fig/{filestem}/upset_{cyto}_context_independent.png"))
upset(fromList(CIG_by_cytokine), order.by = "freq", point.size = 3.5, nintersects = 10,
      line.size = 2, mainbar.y.label = glue("overlap of context \nindependent {cyto} DEGs"),
      sets.x.label = glue("total context \nindependent {cyto} DEGs"),
      text.scale = c(2, 1.7, 1.3, 1.3, 1.7, 2),
      shade.color = getCytoColors()[cyto], sets.bar.color = getCytoColors()[cyto],
      main.bar.color = getCytoColors()[cyto], matrix.color = getCytoColors()[cyto])
#dev.off()

#png(filename = glue("./fig/{filestem}/upset_{cyto}_context_independent.png"))
upset(fromList(CDG_by_cytokine), order.by = "freq", point.size = 3.5, nintersects = 10,
      line.size = 2, mainbar.y.label = glue("overlap of context \ndependent {cyto} DEGs"),
      sets.x.label = glue("total context \ndependent {cyto} DEGs"),
      text.scale = c(2, 1.7, 1.3, 1.3, 1.7, 2),
      shade.color = getCytoColors()[cyto], sets.bar.color = getCytoColors()[cyto],
      main.bar.color = getCytoColors()[cyto], matrix.color = getCytoColors()[cyto])
#dev.off()

### VOLCANO BASED ON CIG/CDG STATUS --------------------------------------------
universal_CIGs <- Reduce(intersect, CIG_by_cytokine)
universal_CDGs <- Reduce(intersect, CDG_by_cytokine)
#CIG_color <- getCytoColors()[cyto]
#CDG_color <- paste0("#",as.hexmode(strtoi(sub("#","",CIG_color),base=16) - 50))
color_mapping <- ifelse(naive_gene_info$X %in% universal_CIGs, "blue",
                        ifelse(naive_gene_info$X %in% universal_CDGs,
                               "red","grey90"))
names(color_mapping)[color_mapping == "blue"] <- "always independent"
names(color_mapping)[color_mapping == "red"] <- "always dependent"
names(color_mapping)[color_mapping == 'grey90'] <- "neither"

EnhancedVolcano(naive_gene_info, lab=NA, x='log2FoldChange', y='padj',
                title=glue("none vs {cyto}"), titleLabSize=16, axisLabSize=12,
                subtitle=NULL, caption=NULL, legendPosition='bottom',
                max.overlaps=20, drawConnectors=TRUE, arrowheads=FALSE,
                colCustom = color_mapping, pCutoff = 0.05)


