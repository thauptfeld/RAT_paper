setwd("C:/Users/acer/OneDrive - Flinders/RATrevision/scripts")
setwd("C:/Users/haup0007/OneDrive - Flinders/RATrevision/scripts")
library(ggplot2)
library(reshape2)
library(Hmisc)
library(dplyr)
library(ape)
library(vegan)


gtdb_nr <- read.csv('../evaluations/20231115.gtdb_vs_nr_gw.csv', header=T)

data<-melt(gtdb_nr)
labels=c('phylum', 'class', 'order', 'family', 'genus', 'species')

ggplot(data = data) + geom_bar(aes(x=variable, y=value, fill=database), 
                               stat='identity', position = position_dodge()) +
  xlab("Rank") + ylab("Fraction of annotated MAGs") +
  scale_x_discrete(labels=labels) +
  theme_bw() 
ggsave('../plots/20231115.gtdb_vs_nr_gw.png', height=4, width=7)


gtdb_nr <- read.csv('../evaluations/20231115.gtdb_vs_nr_cami.csv', header=T)
data<-melt(gtdb_nr)

ggplot(data = data) + geom_bar(aes(x=variable, y=value, fill=database), 
                               stat='identity', position = position_dodge()) +
  xlab("Rank") + ylab("Fraction of annotated MAGs") +
  scale_x_discrete(labels=labels) +
  theme_bw() 
ggsave('../plots/20231115.gtdb_vs_nr_cami.png', height=4, width=7)
