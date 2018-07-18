options(stringsAsFactors = FALSE)

###
library(RColorBrewer)
library(ggplot2)

### Read in data
data <- read.csv("Chemical_reprogramming_trajectory.csv", row.names = 1)

### Set colors to identify cell types
data$Batch <- factor(data$Batch,
                     levels = c("MEF","SI5D","SI12D","XEN","SII8D","SII12D",
                                "SIII3D","SIII6D","SIII8D","SIII10D","SIII15D",
                                "SIII21D"),ordered = TRUE)
myColors <- c("#CC60D8", "#9333B4", "#3701A5", "#120FC0", "#024D63", "#007E00", 
              "#73BB00", "#C3DF23", "#FFF000", "#FED203", "#FF6C01", "#FE0009")
names(myColors) <- c("MEF","SI5D","SI12D","XEN","SII8D","SII12D","SIII3D","SIII6D",
                     "SIII8D","SIII10D","SIII15D","SIII21D")
colScale <- scale_colour_manual(name = c("MEF","SI5D","SI12D","XEN","SII8D","SII12D",
                                         "SIII3D","SIII6D","SIII8D","SIII10D","SIII15D",
                                         "SIII21D"),values = myColors)

### Random the order of cells before plotting
set.seed(1011)
data_plot <- data[sample(rownames(data)),]

### Plot
p <- ggplot() +
  geom_point(data = data_plot, aes(x=data_dim_1, y=data_dim_2, color = Batch), cex = .5) +
  colScale +
  theme_classic() +
  labs(x = "Component 1") +
  labs(y = "Component 2") +
  theme(legend.position="none")

### Output figures
tiff(file=paste0("Chemical_reprogramming_trajectory.tiff"), width = 7, height = 7, units = 'in', 
     res = 300, compression = 'none')
print(p)
dev.off()
