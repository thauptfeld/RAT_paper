setwd("C:/Users/acer/OneDrive - Flinders/RATrevision/scripts")
library(ggplot2)


data <- read.csv('../evaluations/20231028.mousegut_unifrac.csv')

CAMI_levels=c('CAMI.6','CAMI.13','CAMI.23','CAMI.25','CAMI.26',
              'CAMI.30','CAMI.33','CAMI.34','CAMI.38','CAMI.53')
# tool_levels=c('centrifuge','bracken','kraken2','kaiju','robust',
#               'contig','sensitive')
tool_levels=c('centrifuge','bracken','kraken2','kaiju','robust',
              'nobin','wbin','wmetabat')


data$X=factor(data$X,
                     levels=tool_levels)
# 
# axis_labels=c('Centrifuge','Bracken','Kraken2','Kaiju','RAT -mc',
#               'RAT -cr','RAT -mcr')
axis_labels=c('Centrifuge','Bracken','Kraken2','Kaiju','RAT -mc',
              'RAT -cr','RAT -mcr\n(CAMI2)','RAT -mcr\n(MetaBAT2)')

p.data <- ggplot(data, aes(X, variable, fill=value)) + geom_tile() +
  ylab("Sample")+xlab("Reconstructed Profile") + scale_x_discrete(labels=axis_labels) +
  scale_fill_viridis_c(option="D", direction=-1, name='Distance\nWeighted\nUnifrac') + 
  theme(
    axis.text.x=element_text(angle = 270, hjust=0, vjust=0.5),
    text = element_text(size=7),
    legend.box.margin = margin(0,-5,0,-10)
  )

p.data
ggsave('../plots/20231201_Unifrac_heatmap_small_plant.svg', height=75, width=88, units='mm')






