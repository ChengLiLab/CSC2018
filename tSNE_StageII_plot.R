################################################
# Paper can be download from：                 #
# https://www.ncbi.nlm.nih.gov/pubmed/29937202 #
# Cheng Li Lab, Peking University              #
# Qian Zhang(q.zhang.012@gmail.com)            #
################################################
options(stringsAsFactors = F)

### load libraries ----
library(ggplot2)
library(RColorBrewer)

### load data ----
load("Data/data_raw_norm.Robj")
gene_id <- read.csv("Data/fData.csv", stringsAsFactors = FALSE,row.names = 1)
SII_tsne <- read.csv("StageII/AGG_SII_mean0.4_dim10.csv", stringsAsFactors = F, row.names = 1)
set.seed(1011)

### set colors for time points ----
SII_tsne$Batch <- factor(SII_tsne$Batch,
                         levels = c("XEN","SII8D","SII12D"),ordered = TRUE)
myColors <- c("#120FC0", "#024D63", "#007E00")
names(myColors) <- c("XEN","SII8D","SII12D")
colScale <- scale_colour_manual(name = c("XEN","SII8D","SII12D"),values = myColors)

### plot t-SNE with time point information ----
p <- ggplot() +
  geom_point(data = SII_tsne, aes(x=data_dim_1, y=data_dim_2, color = Batch)) +
  colScale +
  theme_classic() +
  # ggtitle("SII") +
  labs(x = " ") + 
  labs(y = " ") + 
  # theme(plot.title = element_text(size = 30, face = "bold", color = "steelblue")) +
  theme(legend.position="none", axis.text=element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

# save figure
tiff(file=paste0("StageII/tSNE_SII_Batch.tiff"), width = 7, height = 7, units = 'in', res = 300, compression = "none")
print(p)
dev.off()

### plot t-SNE with genes expression information ----
data <- SII_tsne[,c(2,3,4,6,10)]

# subset data
barcode <- as.character(data$Barcode)
data_sub <- data_raw[,barcode]
summary(colnames(data_sub) %in% data$Barcode)

### plot t-SNE with gene expression information ----
genes <- c("Sox17", "Sall4", "Oct4", "Zscan4c")
n <- length(barcode)
for (gene in genes) {
  id <- gene_id[gene_id$gene_short_name %in% gene,]
  id <- as.vector(as.matrix(id))
  id <- id[1]
  l <- apply(data_sub[id,] - .1,2,sum) + .1
  f <- l == 0
  l <- log2(l)
  l[f] <- NA
  mi <- min(l,na.rm=TRUE)
  ma <- max(l,na.rm=TRUE)
  data.tsne <- as.data.frame(data.tsne)
  data.tsne$UMI <- l
  l <- is.na(data.tsne$UMI)

  data_grey <- data.tsne[l,]
  data_col <- data.tsne[!l,]
  data_col <- data_col[order(data_col$UMI),]
  data_col <- sample(data_col)
  
# set colors
  myPalette <- colorRampPalette(c(rev(brewer.pal(n = 9,name = "Greys")[1:2]), brewer.pal(n = 9,name = "OrRd"), brewer.pal(n = 9,name = "OrRd")[9]))
  sc <- scale_colour_gradientn(colours = myPalette(256), limits=c(mi, ma))
  p <- ggplot() +
    geom_point(data = data_grey, aes(x=data_dim_1, y=data_dim_2), color = "#F0F0F0", cex = 2.5) +
    geom_point(data = data_col, aes(x=data_dim_1, y=data_dim_2, color = UMI), cex = 2.5) +
    sc +
    theme_void() +
    theme(legend.position=("none"))
  
# save figures for each gene  
  tiff(file=paste0("StageII/SII_tSNE_",gene,"_mean",0.4,"_dim",10,".tiff"), width = 7, height = 7, units = 'in', res = 300, compression = "none")
  print(p)
  dev.off()
}
