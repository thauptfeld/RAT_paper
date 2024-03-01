install.packages('ggalluvial')
library(ggalluvial)

setwd("~/PhD/RAT/scripts")

ucb <- read.csv('../dataset-37099.csv')
ggplot(ucb, aes(y=Freq,axis1=Gender,axis2=Dept)) + geom_alluvium(aes(fill=Admit), width=1/12) +
  geom_stratum(width=1/12, fill='black', color='white') + geom_label(stat='stratum', aes(label=after_stat(stratum)), size=2) +
  scale_x_discrete(limits=c('Gender', 'Dept'), expand=c(.05,.05)) +
  scale_fill_brewer(type='qual', palette='Set1')

ggsave('../plots/test_alluvial.png', width=4, height=4)


rat.new=read.csv('../sankey_new_schem.csv')


ggplot(rat.new, aes(y=freq2_round,axis1=read, axis2=contig, axis3=bin)) +
  geom_alluvium(aes(fill=bin), reverse=F) + geom_stratum(width=1/6, reverse=F) +
  scale_x_continuous(breaks = 1:3, labels = c("Read", 
                                              "Contig", 
                                              "Bin"), expand=c(0,0)) +
  scale_fill_manual(values=c('#A61C3C','#83984d','#f0b51f',
                             '#ea700b','#336b8b','#423f3f')) + 
  ylab('Taxa') +
  theme(
    panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'),
    panel.grid.minor = element_line(colour='#DDDDDD'),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=0.5),
    axis.title.y = element_text(size=20),
    axis.ticks.y = element_blank()
  )
ggsave('../plots/20220902_rat_alluvial_schematic_w_label_new.png', width=6, height=8)





rat.test<- read.csv('../RAT_numbers_schem_smart.csv')

rat.test$Step<-factor(rat.test$Step, levels=c('unclassified','bin', 'contig', 'contig_dm', 'read_dm'))
rat.test$Rank<-factor(rat.test$Rank, levels=c('unclassified','superkingdom', 'phylum', 'class', 
                                              'order','family','genus','species'))
colours=c('#7fb1d3','#5e5959','#ff845e','#f72349','#f6dc85')

ggplot(rat.test, aes(y=newer_freq,axis1=Step,axis2=Rank)) + geom_alluvium(aes(fill=Step, alpha=0.1), alpha=0.8, width=1/12) +
  geom_stratum(width=1/12, fill='black', color='white') + geom_label(stat='stratum', aes(label=after_stat(stratum)), 
                                                                     size=3, label.padding = unit(0.1, 'lines')) +
  scale_x_discrete(limits=c('Step', 'Rank'), expand=c(.07,.07)) +
  scale_fill_manual(values=colours) + theme(
    axis.text.y=element_blank(), axis.ticks.y =element_blank(), panel.background = element_rect(fill='white'),
    panel.grid = element_blank()
  )
ggsave('../plots/rat_alluvial_schematic_w_label.png', width=4, height=4)
