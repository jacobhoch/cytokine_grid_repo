# generate grid of percent DEGs

getPercentDEGs <- function(cyto_A, cyto_B, grid, filestem) {
  
  repol_DEGs <- readtxt(glue("./data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_up.txt"))
  naive_DEGs <- readtxt(glue("./data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_up.txt"))
  
  indep_DEGs <- naive_DEGs[naive_DEGs %in% repol_DEGs]
  percent <- length(indep_DEGs) / length(naive_DEGs)
  
  return(percent)
}

source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-s', '--stem', type='character', nargs=1, help='stem where DEG results live')

args        <- parser$parse_args()
filestem    <- args$s

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
diag(grid) <- 1.0

png(file=glue("./fig/{filestem}_percentDEGgrid.png"))
ht <- Heatmap(as.matrix(grid), 
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", grid[i, j]), x, y, gp = gpar(fontsize = 12))},
              cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE)
draw(ht)
dev.off()

