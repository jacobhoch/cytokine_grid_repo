source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-f', '--filestem', type='character', nargs=1, help='stem where DEG results live')
parser$add_argument('-u', '--uninf_stem', type='character', nargs=1, help='stem where uninf results live')

args          <- parser$parse_args()
filestem      <- args$f
uninf_stem    <- args$u

### IFNB + WT VENN DIAGRAM  ----------------------------------------------------
Mtb_DEGs <- read.csv(file = glue("./data/DE_results/{filestem}/WT_none_vs_none_none_DE.csv"))
Mtb_DEGs <- Mtb_DEGs$X

IFNb_DEGs <- read.csv(file = glue("./data/DE_results/{uninf_stem}/IFNb_none_vs_none_none_DE.csv"))
IFNb_DEGs <- IFNb_DEGs$X

EccCa1_DEGs <- read.csv(file = glue("./data/DE_results/{filestem}/EccCa1.Tn_none_vs_none_none_DE.csv"))
EccCa1_DEGs <- EccCa1_DEGs$X

colors <- c("IFNb"="#CD87F8","Mtb"="#B8D455","EccCa1:Tn"="#31AFF5")
venn_genes <- list("IFNb"=IFNb_DEGs, "Mtb"=Mtb_DEGs, "EccCa1:Tn"=EccCa1_DEGs)

pdf(file = glue("./fig/TBinfection/IFNb_DEG_venn.pdf"),
    height = 5, width = 6)
plot(euler(venn_genes, shape='circle'),
     fill=list(fill=colors[names(venn_genes)], alpha=0.3),
     edges=list(col=colors[names(venn_genes)], lwd=4),
     quantities=list(cex=1.8, col="black", font=2),
     labels=list(cex=2, col=colors[names(venn_genes)]))
dev.off()

### IL4 DEG VENN DIAGRAM --------------------------------===--------------------
MtbIL4_DEGs <- read.csv(file = glue("./data/DE_results/{filestem}/WT_IL4_vs_WT_none_DE.csv"))
MtbIL4_DEGs <- MtbIL4_DEGs$X

IL4_DEGs <- read.csv(file = glue("./data/DE_results/{uninf_stem}/none_IL4_vs_none_none_DE.csv"))
IL4_DEGs <- IL4_DEGs$X

IFNbIL4_DEGs <- read.csv(file = glue("./data/DE_results/{uninf_stem}/IFNb_IL4_vs_IFNb_none_DE.csv"))
IFNbIL4_DEGs <- IFNbIL4_DEGs$X

colors <- c("IFNb\nto IL4"="#CD87F8","Mtb\nto IL4"="#B8D455","IL4"="#929000")
venn_genes <- list("IFNb\nto IL4"=IFNbIL4_DEGs, "Mtb\nto IL4"=MtbIL4_DEGs, "IL4"=IL4_DEGs)

pdf(file = glue("./fig/sfigs/IL4_DEG_venn.pdf"),
    height = 4, width = 6)
plot(euler(venn_genes, shape='circle'),
     fill=list(fill=colors[names(venn_genes)], alpha=0.3),
     edges=list(col=colors[names(venn_genes)], lwd=4),
     quantities=list(cex=1.8, col="black", font=2),
     labels=list(cex=2, col=colors[names(venn_genes)]))
dev.off()


