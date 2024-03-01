setwd("~/PhD/RAT/plots/plot_data")

library(ggplot2)
library(RColorBrewer)
library(svglite)

data=read.csv('wageningen_family_top11.csv')
# tiff("../20230503_wageningen_relabund_top12_family_BIG.tiff", units="in", width=11, height=7, res=600)


ggplot(data, aes(fill=taxon, x=sample, y=value)) + geom_bar(position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired', name='Taxon') + xlab("Sample") + ylab("Relative Abundance") +
  theme(
    axis.text.x = element_text(angle=270, hjust=0.5, vjust=0.5),
    axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill='white'),
    panel.grid = element_blank(), 
    text = element_text(size=7), 
    legend.key.size = unit(4, 'mm')
    # panel.grid.major = element_line(colour="#999999"), panel.grid.minor = element_line(colour='#aaaaaa')
  )

# dev.off()

ggsave('../20230503_wageningen_relabund_top12_family.svg', height=80, width=160, units='mm')
