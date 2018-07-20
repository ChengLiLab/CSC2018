################################################
# Paper can be download from：                 #
# https://www.ncbi.nlm.nih.gov/pubmed/29937202 #
# Cheng Li Lab, Peking University              #
# Qian Zhang(q.zhang.012@gmail.com)            #
################################################
options(stringsAsFactors = F)

### load libraries ----
library(viridis)
library(monocle)

### Load data ----
load("Data/AGG_estimate_size_and_dispersion.Robj")

# Subset Stage II cells
AGG3 <- AGG[, pData(AGG)$Batch %in% c("XEN","SII8D","SII12D")]
AGG3 <- detectGenes(AGG3, min_expr = 0.1)

# Delete Feeders
feeder <- read.csv("Data/feeders in Stage II.csv", stringsAsFactors=FALSE,fileEncoding="utf-16", sep = "\t")
l <- as.character(pData(AGG3)$Barcode) %in% feeder$Barcode
AGG3 <- AGG3[,!l ]

### t-SNE start ----
mean <- 0.4
disp_table <- dispersionTable(AGG3)
ordering_genes <- subset(disp_table,
                         mean_expression >= mean &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id
n <- length(ordering_genes)
AGG3 <- setOrderingFilter(AGG3, ordering_genes)
plot_ordering_genes(AGG3)

dim <- 10
AGG3 <- AGG3
AGG3 <- reduceDimension(AGG3, max_components = 2, norm_method = 'log',
                            reduction_method = 't-SNE', num_dim = dim, verbose = T) # Change 2 num_dim 2..15
AGG3 <- clusterCells(AGG3, num_clusters = 5, verbose = T)

### save t-SNE dimensions ----
cds <- AGG3
x=1
y=2
color_by="Cluster"
markers=NULL
show_cell_names=FALSE
cell_size=1.5
cell_name_size=2
if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
  stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
}

gene_short_name <- NULL
sample_name <- NULL
data_dim_1 <- NULL
data_dim_2 <- NULL

lib_info <- pData(cds)

t-SNE_dim_coords <- reducedDimA(cds)
data_df <- data.frame(t(t-SNE_dim_coords[c(x,y),]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- colnames(cds)
data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")

# save dimensions
write.csv(data_df, file = paste0("StageII/AGG_SII_mean",mean,"_dim",dim,".csv"), quote = F)
