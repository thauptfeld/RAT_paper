library(ggplot2)
library(scatterpie)

setwd("C:/Users/acer/OneDrive - Flinders/RATrevision/scripts")
levels=c('diamond', 'r_nobin', 'r_wbin', 'r_wmetabat', 'r_robust')

mouses_scatter=read.csv('../evaluations/20231027.CAMI2_mousegut_scatterpie_onlyRAT.csv', header=T)
mouses_scatter$tool<-factor(mouses_scatter$tool, levels=levels)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.05), data=mouses_scatter,
                           cols=c("avg_classified", "avg_unclassified"),size=0.3) + 
  theme_bw() +
  scale_x_continuous(breaks = c(1.3, 1.9, 2.5), labels=c('phylum', 'genus', 'species')) +
  xlab("Taxonomic rank") + ylab("Mean TPR") + coord_fixed() +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(text=element_text(size=7),
        legend.key.size = unit(4, 'mm'))

ggsave('20240227.correct_incorrect_mousegut_onlyRAT.svg', height=80, width=88, units = 'mm')


marine_scatter=read.csv('../evaluations/20231024.CAMI2_marine_scatterpie_final.csv', header=T)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.05), data=marine_scatter,
                           cols=c("avg_A", "avg_B"),size=0.3) + 
  theme_bw() +
  scale_x_continuous(breaks = c(1.4, 2.2, 3), labels=c('phylum', 'genus', 'species')) +
  xlab("Taxonomic rank") + ylab("Average TPR") + coord_fixed() +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('../plots/20231027.correct_incorrect_marine_scatterpie_draft.png', height=4, width=5)


plant_scatter=read.csv('../evaluations/20231024.CAMI2_plant_scatterpie_3_ranks.csv', header=T)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.05), data=plant_scatter,
                           cols=c("avg_A", "avg_B"),size=0.3) + 
  theme_bw() +
  scale_x_continuous(breaks = c(0.4, 1.2, 2), labels=c('phylum', 'genus', 'species')) +
  xlab("Taxonomic rank") + ylab("Average TPR") + coord_fixed() +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('../plots/20231027.correct_incorrect_plant_scatterpie_plant_draft.png', height=4, width=5)

mouse_scatter=read.csv('../evaluations/20231027.CAMI2_mousegut_scatterpie_3_ranks.csv', header=T)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.05), data=mouse_scatter,
                           cols=c("avg_A", "avg_B"),size=0.3) + 
  theme_bw() +
  scale_x_continuous(breaks = c(1.4, 2.3, 3.2), labels=c('phylum', 'genus', 'species')) +
  xlab("Taxonomic rank") + ylab("Mean TPR") + coord_fixed() +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('../plots/20231027.correct_incorrect_scatterpie_mousegut_draft.png', height=4, width=5)


all_small=rbind(marine_scatter, plant_scatter, mouses_scatter)
all_small$environment <- gsub("plant", "rhizosphere", all_small$environment)
all_small$environment=factor(all_small$environment, levels=c("mouse gut", "marine", "rhizosphere"))



ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.05), data=all_small,
                           cols=c("avg_A", "avg_B"),size=0.3) + 
  theme_bw() +
  scale_x_continuous(breaks = c(1.45,2.35,3.25), labels=c("phylum","genus","species")) +
  xlab("Taxonomic rank") + ylab("Average TPR") + coord_fixed() +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(legend.position = "none") + facet_wrap(~environment, ncol=1)
ggsave('../plots/20231031.correct_incorrect_scatterpie_allenv_3_ranks_draft.png', height=7, width=5)


br=c(0.45,1.35,2.25,3.15,4.05,4.95)
l=c('phylum', "class", "order", "family", 'genus', 'species')



all_scatter=read.csv('../evaluations/20231027.CAMI2_allenv_scatterpie_all_ranks.csv', header=T)

all_scatter$environment=factor(all_scatter$environment, levels=c("mouse gut", "marine", "plant"))

ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.05), data=all_scatter,
                           cols=c("avg_A", "avg_B"),size=0.3) + 
  theme_bw() +
  scale_x_continuous(breaks = br, labels=l) +
  xlab("Taxonomic rank") + ylab("Mean TPR") + coord_fixed() +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(legend.position = "none") + facet_wrap(~environment, ncol=1)
ggsave('../plots/20240227.correct_incorrect_scatterpie_allenv_ar_draft.svg', height=150, width=185, units = 'mm')
