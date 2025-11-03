# generate buckets of DEGs
# then maybe make a heatmap out of them? who's to say

getDEcategory <- function(naive_genes, repol_DEGs) {
  row <- c()
  row[1] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE", 1])
  row[2] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE likely", 1])
  row[3] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE not excluded", 1])
  row[4] <- sum(naive_genes %in% repol_DEGs[repol_DEGs$category == "DE unlikely", 1])
  #print(row)
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

repol_DEGs <- read.csv(glue("./data/DE_results/{filestem}/{cyto_A}_{cyto_B}_vs_{cyto_A}_none_full.csv"))
naive_DEGs <- read.csv(glue("./data/DE_results/{filestem}/none_{cyto_B}_vs_none_none_full.csv"))

DE_both <- repol_DEGs$X[repol_DEGs$DE & !repol_DEGs$notDE & naive_DEGs$DE & !naive_DEGs$notDE]
notDE_then_DE <- repol_DEGs$X[repol_DEGs$DE & !repol_DEGs$notDE & !naive_DEGs$DE & naive_DEGs$notDE]
DE_then_notDE <- repol_DEGs$X[!repol_DEGs$DE & repol_DEGs$notDE & naive_DEGs$DE & !naive_DEGs$notDE]
DE_then_cusp <- repol_DEGs$X[!repol_DEGs$DE & !repol_DEGs$notDE & naive_DEGs$DE & !naive_DEGs$notDE]

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
#plot(tableGrob(grid, theme=ttheme_default(base_size=25)))
print(grid)
