library(ggplot2)
library(scatterpie)

setwd("C:/Users/acer/OneDrive - Flinders/RATrevision/scripts")

marine_scatter=read.csv('../revision/evaluations/20231024.CAMI2_marine_scatterpie_final.csv', header=T)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.03), data=marine_scatter,
                           cols=c("avg_A", "avg_B"),size=0.2) + theme_bw() +
  scale_x_continuous(breaks = c(1.4, 2.3, 3.2), labels=c('phylum', 'genus', 'species'))
ggsave('20231024.correct_incorrect_scatterpie_draft.png', height=2.5, width=7)


plant_scatter=read.csv('../revision/evaluations/20231024.CAMI2_plant_scatterpie_final.csv', header=T)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.03), data=plant_scatter,
                           cols=c("avg_A", "avg_B"),size=0.2) + theme_bw() +
  scale_x_continuous(breaks = c(1.4, 2.3, 3.2), labels=c('phylum', 'genus', 'species'))
ggsave('20231024.correct_incorrect_scatterpie_plant_draft.png', height=2.5, width=7)

mouses_scatter=read.csv('../evaluations/20231027.CAMI2_mousegut_scatterpie_onlyRAT.csv', header=T)
ggplot() + geom_scatterpie(aes(x=x.axis, y=avg_TPR, r=0.03), data=mouses_scatter,
                           cols=c("avg_A", "avg_B"),size=0.2) + theme_bw() +
  scale_x_continuous(breaks = c(1.4, 2.3, 3.2), labels=c('phylum', 'genus', 'species'))
ggsave('20231024.correct_incorrect_scatterpie_plant_draft.png', height=2.5, width=7)
