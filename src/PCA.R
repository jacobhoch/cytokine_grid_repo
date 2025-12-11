source("./src/utilities.R")
prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-d', '--dds_list', type='character', nargs=1, help='list of DESeq objects')
parser$add_argument('-c', '--components', type='character', nargs='*', default='1,2', help='')
parser$add_argument('-f', '--filestem', type='character', nargs=1, help='')

args        <- parser$parse_args()
dds_list    <- args$d
components  <- args$c
filestem    <- args$f

cytokines <- getCytokines()
cyto_colors <- getCytoColors()

set.seed(123)
dds <- dds_list[[1]]

# borrowed from DESeq source code, but I want to save intermediates
vsd <- vst(dds, blind=FALSE, fitType = 'local')
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Donor)

rv <- rowVars(assay(vsd)) # calculate the variance for each gene
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd)[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
# add the intgroup factors together to create a new grouping factor
intgroup <- c("Stim","day_1","day_2","Donor","tech_rep")
intgroup <- c("Donor")
group <- if (length(intgroup) > 1) {
  intgroup.df <- as.data.frame(colData(vsd)[, intgroup, drop=FALSE])
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(vsd)[[intgroup]]
}

# assembly the data for the plot
#barplot(percentVar[1:10], xlab="component", ylab="percent variance")
pcsToUse <- 1:3
pcs <- paste0("PC", pcsToUse)
pcaData <- data.frame(V1=pca$x[,pcsToUse[1]],
                      V2=pca$x[,pcsToUse[2]],
                      V3=pca$x[,pcsToUse[3]],
                      group=group, name=colnames(vsd), colData(vsd))
colnames(pcaData)[1:3] <- pcs
attr(pcaData, "percentVar") <- percentVar[pcsToUse]

percentVar <- round(100 * attr(pcaData, "percentVar"))

# combine day_1 and day_2 splits into one dataframe
pca_metadata <- cbind(as.character(vsd$day_1),as.character(vsd$day_2))
rownames(pca_metadata) <- colnames(vsd)
pcaData <- merge(pcaData, pca_metadata, by.x = "name", by.y = "row.names")

cols <- c()
for (i in cytokines) {
  pcaData[[i]] <- with(pcaData, ifelse(day_2 == i, 1, 0))
  cols <- c(cols,i)
}

pcaData <- pivot_longer(pcaData,cols=all_of(cols),
                        names_to="colorby",values_to='amount')

png(filename = glue(".fig/{filestem}_PCA.png"))
ggplot(pcaData) +
  geom_arc_bar(aes(x0=PC1,y0=PC2,r0=0.3,r=1,amount=amount,
                   fill=colorby,col=colorby), stat='pie', linewidth=0) +
  coord_equal() +
  theme_classic() +
  scale_fill_manual(name='cytokine 2',values=cyto_colors) +
  geom_point(aes(x=PC1,y=PC2,color=day_1), size=4.5) +
  scale_color_manual(name='cytokine 1',values=colors) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  guides(fill = guide_legend(override.aes = list(linetype = rep(1,7)))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()
