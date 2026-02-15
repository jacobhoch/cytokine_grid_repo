# generate grid of percent DEGs

# getPercentDEGs <- function(cyto_A, cyto_B, grid, filestem) {
#   
#   repol_DEGs <- readtxt(glue("./data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_up.txt"))
#   naive_DEGs <- readtxt(glue("./data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_up.txt"))
#   
#   indep_DEGs <- naive_DEGs[naive_DEGs %in% repol_DEGs]
#   percent <- length(indep_DEGs) / length(naive_DEGs)
#   
#   return(percent)
# }

getPercentDEGs <- function(cyto_A, cyto_B, grid, filestem) {
  
  repol_DEGs_up <- readtxt(glue("./data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_up.txt"))
  naive_DEGs_up <- readtxt(glue("./data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_up.txt"))
  
  repol_DEGs_down <- readtxt(glue("./data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_down.txt"))
  naive_DEGs_down <- readtxt(glue("./data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_down.txt"))
  
  indep_DEGs <- naive_DEGs_up[naive_DEGs_up %in% repol_DEGs_up]
  indep_DEGs <- c(indep_DEGs, naive_DEGs_down[naive_DEGs_down %in% repol_DEGs_down])
  percent <- length(indep_DEGs) / length(c(naive_DEGs_up, naive_DEGs_down))
  
  return(percent)
}

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-f', '--filestem', type='character', nargs=1, help='stem where DEG results live')

args        <- parser$parse_args()
filestem    <- args$f

cytokines <- getCytokines()

grid <- data.frame(matrix(nrow = length(cytokines), ncol = length(cytokines)),
                   row.names = cytokines)
colnames(grid) <- cytokines

conditions <- unique(str_split_i(list.files(glue("./data/DE_results/{filestem}")), pattern='_vs_', i=1))
conditions <- conditions[!grepl("none",conditions)]

for (i in 1:length(conditions)){
  # pull cytokine information to determine which DEGs to look at
  stim <- unlist(str_split(conditions[i], pattern='_'))
  cyto_A <- stim[1]
  cyto_B <- stim[2]
  
  grid[cyto_A,cyto_B] <- getPercentDEGs(cyto_A, cyto_B, grid, filestem)
}
diag(grid) <- 1
colors <- colorRampPalette(brewer.pal(9, "BuPu"))(255)

pdf(file=glue("./fig/cytokine_grid/indep_DEG_grid.pdf"),
    width = 7, height = 5)
ht <- Heatmap(as.matrix(grid), col=colors, rect_gp = gpar(col = "black", lwd=1.5),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", grid[i, j]), x, y, gp = gpar(fontsize = 16, col="black", fontface = "bold"))},
              cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = TRUE,
              row_title = "context (day 1)", column_title = "cytokine (day 2)",
              column_title_side = "bottom", width = unit(10, "cm"), height = unit(10, "cm"),
              row_names_side = "left", na_col = "black",
              heatmap_legend_param = list(title="context dependence\nscore",
                                          title_gp=gpar(fontsize=12),
                                          labels_gp=gpar(fontsize=10),
                                          legend_height=unit(3,'cm'),
                                          legend_width=unit(1,'cm')))
draw(ht)
dev.off()

