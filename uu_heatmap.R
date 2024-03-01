setwd("~/PhD/RAT")
library(ggplot2)

data <- read.csv('distance_matrix_melted_only_important.csv')
data <- read.csv('distance_matrix_melted_no_CAMI_new.csv')

CAMI_levels=c('CAMI.6','CAMI.13','CAMI.23','CAMI.25','CAMI.26',
              'CAMI.30','CAMI.33','CAMI.34','CAMI.38','CAMI.53')
tool_levels=c('centrifuge','bracken','kraken2','kaiju','robust',
              'nobin','wbin','wmetabat')

data$variable=factor(data$variable, 
                        levels=rev(CAMI_levels))
data$X=factor(data$X, 
                     levels=tool_levels)

axis_labels=c('Centrifuge','Bracken','Kraken2','Kaiju','RAT robust',
              'RAT without MAGs','RAT CAMI genomes','RAT MetaBAT2 MAGs')

p.data <- ggplot(data, aes(X, variable, fill=value)) + geom_tile() +
  ylab("True Profile")+xlab("Reconstructed Profile") + scale_x_discrete(labels=axis_labels) +
  scale_fill_viridis_c(option="D", direction=-1, name='Distance\nWeighted\nUnifrac') + 
  theme(
    axis.text.x=element_text(angle = 270, hjust=0, vjust=0.5),
    plot.title = element_text(size=20)
  )
p.data
ggsave('plots/20230227_Unifrac_heatmap_small_weighted.png', height=4.5, width=4.5)



setwd("~/PhD/RAT/unifrac")
library(ggplot2)
library(gplots)

uu_input<-read.csv('distance_matrix_small_correct.csv', row.names=1)
uu_dist<-dist(uu_input)
uu_clust<-hclust(uu_dist)
plot(uu_clust)

mycolors <- colorRampPalette(c("blue", "green", "yellow"))
heatmap.2(t(as.matrix(uu_input), trace = "none", col = mycolors, Rowv = F)


