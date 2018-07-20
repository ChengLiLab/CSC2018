################################################
# Paper can be download from：                 #
# https://www.ncbi.nlm.nih.gov/pubmed/29937202 #
# Cheng Li Lab, Peking University              #
# Qian Zhang(q.zhang.012@gmail.com)            #
################################################
options(stringsAsFactors = F)

### load libraries -----
library(monocle)

### load data ----
# Due the limitation of Github, you should first download
# processed data GSE114952_AGG.tar.gz from GEO:GSE114952.
load("Data/matrix.RData")
pd <- read.csv("Data/pData.csv", row.names = 1)
fd <- read.csv("Data/fData.csv", row.names = 1)

### pre-processing ----
AGG <- newCellDataSet(matrix,
                      phenoData = new("AnnotatedDataFrame", pd),
                      featureData = new("AnnotatedDataFrame", fd),
                      lowerDetectionLimit=0.5,
                      expressionFamily=negbinomial.size())

# estimate_size_and_dispersion, eval=TRUE
AGG <- estimateSizeFactors(AGG)
AGG <- estimateDispersions(AGG, cores = 3)
AGG <- detectGenes(AGG, min_expr = 0.1)

# export raw data
data_raw <- exprs(AGG)
data_raw <- as.matrix(data_raw)
dim(data_raw)

# normalize data
data_norm <- sweep(data_raw,2,colSums(data_raw),'/')*median(colSums(data_raw))

# Convert to data frame
data_raw <- as.data.frame(data_raw)
data_norm <- as.data.frame(data_norm)

### save results ----
save(AGG, file = "Data/AGG_estimate_size_and_dispersion.Robj")
save(data_raw, data_norm, file = "Data/data_raw_norm.Robj")
sessionInfo()
