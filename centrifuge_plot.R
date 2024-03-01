setwd("C:/Users/haup0007/OneDrive - Flinders/RATrevision/checkdata")

library(ggplot2)
library(ggpubr)

data<-read.csv('centrifuge_values_normalized.csv')


colors=c("#66ccee", "#CCCCCC", "#888888", "#333333")
labels=c('Groundwater*', 'Marine', 'Mouse gut', 'Rhizosphere')

line<- ggplot(data) + geom_line(aes(x=matchLength, y=counts, color=environment), size=1) +
  xlim(0,300) + scale_y_continuous(trans='log10') +
  xlab("Database Match Length") + ylab("Fraction of reads (log10)") +
  scale_color_manual(values=colors, labels=labels, name='Environment') +
  theme_bw() 

ggsave("../plots/20231122_centrifuge_matchLengths_normalized.png", width=6.5, height=4)


line<-ggplot(data) + geom_line(aes(x=matchLength, y=counts, color=environment), size=1.5) +
  xlim(0,300) + scale_y_continuous(trans='log10') +
  xlab("Database Match Length") + ylab("Fraction of reads (log10)") +
  scale_color_manual(values=colors, labels=labels) +
  theme_bw() + theme(text=element_text(size=25), legend.position = c(0.84, 0.84),
                     legend.title = element_blank(), legend.key.size = unit(1.3,'cm'),
                     legend.box.background = element_rect(color = "black"))

ggsave("../plots/20231122_centrifuge_matchLengths_normalized_BIG.png", width=11, height=9)



setwd("C:/Users/haup0007/OneDrive - Flinders/RATrevision/")

library(ggplot2)

data<-read.csv('simulated_data_barplot.csv')


colors=c("#227788", "#66ccee", "#fc9843", "#5a1800")
colors=c("#66ccee", "#CCCCCC", "#888888", "#333333")
labels=c('Groundwater*', 'Marine', 'Mouse gut', 'Rhizosphere')

mags<-ggplot(data, aes(x=environment, y=mean_bin, fill=environment, width=1)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Environment") + ylab("Fraction of reads in MAGs") +
  scale_x_discrete(labels=labels) +
  geom_errorbar(aes(ymin=mean_bin-sd_bin, ymax=mean_bin+sd_bin), width=0.2,
                ) +
  scale_fill_manual(values=colors, labels=labels, name='Environment') +
  theme_bw() + theme(text=element_text(size=25), legend.position="none",
                     axis.text.x = element_blank())


ggsave("plots/20231123_environments_MAG_reads_BIG.png", width=4, height=5)

colors=c("#CCCCCC", "#888888", "#333333")
data2<-data[c(1:3),]
unknown<-ggplot(data2, aes(x=environment, y=mean_genus, fill=environment, width=1)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Environment") + ylab("Unknown species") +
  scale_x_discrete(labels=labels) +
  geom_errorbar(aes(ymin=mean_genus-sd_genus, ymax=mean_genus+sd_genus), width=0.2,
  ) +
  scale_fill_manual(values=colors, labels=labels, name='Environment') +
  theme_bw() + theme(text=element_text(size=25), legend.position = "none",
                     axis.text.x = element_blank())


ggsave("plots/20231123_environments_genus_reads_BIG.png", width=4, height=5)


figure <- ggarrange(line,
                    ggarrange(mags, unknown,
                              ncol = 1, nrow = 2), 
                    ncol=2, widths=c(3,1))
ggexport(figure, filename='plots/20231130.suppfig2.png', height=900, width=1600)
