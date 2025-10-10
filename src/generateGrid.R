# as a lapply:
# list of all stims
# function that takes one condition: "A" to "B"
# pulls DEGs A to B vs A to none
# pulls DEGs none to B vs none to none
# report percentage

getPercentDEGs <- function(cyto_A, cyto_B, grid, filestem) {
  
  repol_DEGs <- readtxt(glue("./data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_up.txt"))
  naive_DEGs <- readtxt(glue("./data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_up.txt"))
  
  indep_DEGs <- naive_DEGs[naive_DEGs %in% repol_DEGs]
  percent <- length(indep_DEGs) / length(naive_DEGs)
  
  return(percent)
}

source("./src/utilities.R")
prepEnv()

cytokines <- getCytokines()

grid <- data.frame(matrix(nrow = length(cytokines), ncol = length(cytokines)),
                   row.names = cytokines)
colnames(grid) <- cytokines

test <- unique(str_split_i(list.files(glue("./data/DE_results/{filestem}")), pattern='_vs_', i=1))
test <- test[!grepl("none",test)]

for (i in 1:length(test)){
  # pull cytokine information to determine which DEGs to look at
  stim <- unlist(str_split(test[i], pattern='_'))
  cyto_A <- stim[1]
  cyto_B <- stim[2]
  
  grid[cyto_A,cyto_B] <- getPercentDEGs(cyto_A, cyto_B, grid, filestem)
}
diag(grid) <- 1.0
Heatmap(as.matrix(grid), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", grid[i, j]), x, y, gp = gpar(fontsize = 12))},
        cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE)
