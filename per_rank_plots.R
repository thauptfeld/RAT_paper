library(ggplot2)

setwd("~/PhD/RAT")
setwd("C:/Users/haup0007/OneDrive - Flinders/RATrevision/scripts")

data<-read.csv('../reads_per_step/20231201_per_rank_classified_gtdb_nr.csv')

data$rank <- factor(data$rank, levels=c("superkingdom", "phylum", "class", "order",
                                            "family", "genus", "species"))
data$database <- factor(data$database, levels=c("gtdb", "nr"))

ggplot(data, aes(x=rank, y=classified, fill=database)) +
  geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=c('#f8766d','#00bfc4'), name="Database", labels=c("GTDB", "nr")) +
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major = element_line(colour='#BBBBBB'),
        panel.grid.minor = element_line(colour='#DDDDDD'),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) + 
  xlab('Taxonomic Rank') + ylab('Fraction of annotated reads')
  
ggsave('../plots/20231201_per_rank_classified_gtdb_nr.png', height=4, width=6)
