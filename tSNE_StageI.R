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

### Subset SI12D and XEN cells ----
AGG2 <- AGG[, pData(AGG)$Batch %in% c("SI12D","XEN")] 
AGG2 <- detectGenes(AGG2, min_expr = 0.1)

### t-SNE start ----
mean <- 0.5
disp_table <- dispersionTable(AGG2)
ordering_genes <- subset(disp_table,
                         mean_expression >= mean &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id

AGG2 <- setOrderingFilter(AGG2, ordering_genes)
plot_ordering_genes(AGG2)

dim <- 7
AGG2 <- AGG2
AGG2 <- reduceDimension(AGG2, max_components = 2, norm_method = 'log',
                            reduction_method = 'tSNE', num_dim = dim, verbose = T) # Change 2 num_dim 2..15
AGG2 <- clusterCells(AGG2, num_clusters = 3, verbose = T)

### save t-SNE dimensions ----
cds <- AGG2
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

tSNE_dim_coords <- reducedDimA(cds)
data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- colnames(cds)
data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")

# save dimensions
write.csv(data_df, file = paste0("StageI/AGG_SI12D_XEN_mean",mean,"_dim",dim,".csv"), quote = F)