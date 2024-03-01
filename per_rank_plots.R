library(ggplot2)

setwd("~/PhD/RAT")

data<-read.csv('per_rank_classified.csv')

data$rank <- factor(data$rank, levels=c("superkingdom", "phylum", "class", "order",
                                            "family", "genus", "species"))
data$Category <- factor(data$Category, levels=c("simulated", "biological"))

ggplot(data, aes(x=rank, y=value, fill=Category)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major = element_line(colour='#BBBBBB'),
        panel.grid.minor = element_line(colour='#DDDDDD'),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) + 
  xlab('Taxonomic Rank') + ylab('Fraction of classified reads')
  
ggsave('plots/per_rank_classified_cami_gw.png', height=4, width=6)
