# upset plot of all new DEG categories given previous category
source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')
parser$add_argument('-c', '--cyto', type='character', nargs=1, help='main cytokine')

args        <- parser$parse_args()
filestem    <- args$s
cyto        <- args$c

naive_gene_info <- read.csv(glue("./data/DE_results/{filestem}/none_{cyto}_vs_none_none_full.csv"))

naive_DEGs <- naive_gene_info[naive_gene_info$category == "DE",]$X

cytokines <- getCytokines()

CDG_by_cytokine <- list()
CIG_by_cytokine <- list()
for (i in cytokines[cytokines != cyto]) {
  cyto_gene_info <- read.csv(glue("./data/DE_results/{filestem}/{i}_{cyto}_vs_{i}_none_full.csv"))
  
  CIGs <- intersect(naive_DEGs, cyto_gene_info[cyto_gene_info$category == "DE",]$X)
  CDGs <- intersect(naive_DEGs, cyto_gene_info[cyto_gene_info$category != "DE",]$X)
  
  CIG_by_cytokine[[i]] <- CIGs
  CDG_by_cytokine[[i]] <- CDGs
}

pdf(file = glue("./fig/cytokine_grid/upset_{cyto}_ind.pdf"),
    width = 6, height = 9, onefile = FALSE)
upset(fromList(CIG_by_cytokine), order.by = "freq", point.size = 3.5, nintersects = 10,
      line.size = 2, mainbar.y.label = glue("overlap of context \nindependent {cyto} DEGs"),
      sets.x.label = glue("total context \nindependent {cyto} DEGs"),
      text.scale = c(2, 1.7, 1.3, 1.3, 1.7, 2),
      shade.color = getCytoColors()[cyto], sets.bar.color = getCytoColors()[cyto],
      main.bar.color = getCytoColors()[cyto], matrix.color = getCytoColors()[cyto])
dev.off()

pdf(file = glue("./fig/cytokine_grid/upset_{cyto}_dep.pdf"),
    width = 6, height = 9, onefile = FALSE)
upset(fromList(CDG_by_cytokine), order.by = "freq", point.size = 3.5, nintersects = 10,
      line.size = 2, mainbar.y.label = glue("overlap of context \ndependent {cyto} DEGs"),
      sets.x.label = glue("total context \ndependent {cyto} DEGs"),
      text.scale = c(2, 1.7, 1.3, 1.3, 1.7, 2),
      shade.color = getCytoColors()[cyto], sets.bar.color = getCytoColors()[cyto],
      main.bar.color = getCytoColors()[cyto], matrix.color = getCytoColors()[cyto])
dev.off()

