library(ggplot2)
library(reshape2)
require(gridExtra)

setwd('~/PhD/RAT/scripts')

data.16S <- read.csv('../wageningen_16S_distance_matrix_for_hclust.csv', row.names=1)
data.mgx <- read.csv('../wageningen_distance_matrix_0.001_for_hclust.csv', row.names=1)

test.16S<-as.matrix(data.16S)
test.mgx<-as.matrix(data.mgx[1:17, 1:17])

png("../dendro_heatmap_16s.png")
heatmap.2(test.16S, trace="none", scale="none")
dev.off()
png("../dendro_heatmap_mgx.png")
heatmap.2(test.mgx, trace="none", scale="none")
dev.off()

melt.16S <- reshape2::melt(data.16S)
melt.mgx <- reshape2::melt(data.mgx[1:17,1:18])

install.packages("gplots")
library(gplots)

png("../plots/16S_clustering_avg.png")
plot(hclust(as.dist(test.16S), method="average"))
dev.off()

p.16S <- ggplot(melt.16S, aes(X, variable, fill=value)) + geom_tile() +
  ylab("")+xlab("") + ggtitle("16S Amplicon Sequences") +
  scale_fill_viridis_c(option="D", direction=-1) + 
  theme(
    axis.text.x=element_text(angle = 270), plot.title = element_text(size=20)
  )
p.mgx <- ggplot(melt.mgx, aes(X, variable, fill=value)) + geom_tile() +
  ylab("")+xlab("") + ggtitle("WGS Metagenomics + RAT") +
  scale_fill_viridis_c(option="D", direction=-1, limits = c(0.0, 0.4), oob = scales::squish)  + 
  theme(
    axis.text.x=element_text(angle = 270), plot.title = element_text(size=20)
  )

grid.arrange(p.16S, p.mgx, ncol=2)

ggsave("../plots/wunifrac_comparison_0.1.png", arrangeGrob(p.16S, p.mgx, ncol=2), width=14)
