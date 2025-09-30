format_fc <- function(rds_filename, barcodes_filename) {
  
  AllCounts <- readRDS(rds_filename)
  dge <- as.matrix(AllCounts$umicount$exon$all) # pulling out exon counts
  barcodes_file <- read.csv(barcodes_filename)
  metadata <- data.frame(barcodes_file)
  
  fc <- c()
  fc.colnames <- c()
  for (i in 1:length(barcodes_file$BC)) {
    fc <- cbind(fc,dge[,barcodes_file$BC[i] == colnames(dge)])
    if (barcodes_file$BC[i] %in% colnames(dge)) {
      fc.colnames <- c(fc.colnames, barcodes_file$Unique[i])
    } else {
      print(paste0("missing ",barcodes_file$Unique[i]," (",
                   barcodes_file$Name[i],")"))
      metadata <- metadata[-i,]
    }
  }
  colnames(fc) <- fc.colnames
  
  return(list(as.data.frame(fc), metadata))
}

source('src/utilities.R')


prepEnv()

parser <- ArgumentParser(description="file i/o for oading counts data")
parser$add_argument('-i', '--input_file', type='character', nargs='+', help='input counts data')
parser$add_argument('-m', '--metadata', type='character', nargs='+', help='sample metadata') 
parser$add_argument('-s', '--stem', type='character', nargs='+', help='file stem for saving')
parser$add_argument('-d', '--design', type='character', nargs='+', help='design variables to combine, seperated by spaces')

args        <- parser$parse_args()
input_files <- args$i
sample_info <- args$m
design      <- args$d
stem        <- args$s

# loading gene name stuff
hub = AnnotationHub()
query(hub, c("EnsDb", "Homo sapiens", "97"))
edb = hub[["AH73881"]]
keys = keys(edb, "GENENAME")
columns =  c("GENEID", "ENTREZID", "GENEBIOTYPE")
tbl =
  ensembldb::select(edb, keys, columns, keytype = "GENENAME") %>%
  as_tibble()

# making feature counts matrices for every sequencing run/library
fc <- data.frame()
metadata <- data.frame()
for (i in 1:length(input_files)) {
  barcodes_file <- sample_info[i]
  rds_file <- input_files[i]
  
  c(fc.temp, metadata.temp) %<-% format_fc(rds_file, barcodes_file)
  
  # cleaning up fc and merging into one matrix
  fc <- merge(fc, fc.temp, by.x='row.names', by.y='row.names', all=TRUE)
  rownames(fc) <- fc$Row.names
  fc <- fc[,-1]
  
  # merging the metadata too
  metadata <- rbind(metadata, metadata.temp)
}

#replace GENEID with GENENAME
fc <- cbind(sub('\\.[0-9]*$','', rownames(fc)), fc)
colnames(fc)[1] <- "geneid"
fc <- merge(fc, tbl[,1:2], by.x="geneid", by.y="GENEID", all.y =TRUE)

# reformat for DEseq
fc <- fc[, c("GENENAME",
             names(fc)[names(fc) != "GENENAME"])]
fc <- fc[,-2]
colnames(fc)[1] <- 'gene'
fc[,-1] <- lapply(fc[,-1], function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})
fc[,1] <- make.names(fc[,1], unique=TRUE)
#high_count_genes <- fc$gene[rowSums(fc[,-1], na.rm=TRUE) >= 10]
#countData <- fc[fc$gene %in% high_count_genes,]
#countData[is.na(countData[,])] <- 0

fc[is.na(fc[,])] <- 0
fc <- fc[rowSums(fc[-1]) > 0, colSums(fc[-1]) > 0]
rownames(metadata) <- colnames(fc)[-1]
metadata <- data.frame(lapply(metadata,factor))

design <- as.formula(paste0("~",paste(design, collapse="+")))
dds <- DESeqDataSetFromMatrix(countData=fc,
                              colData=metadata,
                              design=design, tidy=TRUE)

outstem = './data/raw_dds/'
outext  = '.Rds'

saveRDS(dds, glue('{outstem}{stem}{outext}'))