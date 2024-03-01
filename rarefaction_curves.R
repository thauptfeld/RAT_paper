setwd("~/PhD/RAT/rarefaction/CAMI")
library(ggplot2)
library(dplyr)
library(svglite)

smp.ckrr<-read.csv('../smp13.ckrr.wmin.nomin.csv')
tool_colors=c('#b9f291','#f6e5ad','#5e5959','#f72349')
smp.ckrr$tool <- factor(smp.ckrr$tool, levels=c("kaiju", "rat_dm", 'rat_rob',"cami"))
sizes=c(2,2,2,1)

p.rat<-ggplot(smp.ckrr, aes(x=percent, y=taxa, colour=tool, shape=min)) + geom_point(size=3) + 
  stat_summary(geom = "line", fun = mean, size=2) + #geom_hline(yintercept =128, color='#f72349', size=0.1) +
  scale_y_continuous(trans='log10') + 
  scale_colour_manual(values=tool_colors, name='Tool',
                      labels=c('Kaiju', 'RAT -mcr', 'RAT -mc', 'CAMI reference')) + 
  scale_shape_discrete(labels=c('All taxa', 'Taxa>=0.001%'), name='') +
  scale_size_manual(values=sizes) +
  xlab('Percent reads sampled') + ylab('Number of detected taxa') +
  theme(
    panel.background = element_rect(fill='white', colour='#d9d9d9'),
    panel.grid.major = element_line(colour='#bbbbbb'),
    panel.grid.minor = element_line(colour='#dddddd'),
    legend.key.size = unit(1.5,"cm"),
    text = element_text(size=30)
  )

p.rat
ggsave('../../plots/rarefaction/20231031.smp13.ckrr.wmin.nomin.png', width=12, height=7)







smp.krr<-read.csv('../W23-2.krr.wmin.nomin.csv')
tool_colors=c('#b9f291','#f6e5ad','#5e5959')
smp.krr$tool <- factor(smp.krr$tool, levels=c("kaiju", "rat", 'rat_r'))

#tiff("../../plots/rarefaction/20231031.W23-2.krr.wmin.nomin.tiff", units="in", width=9, height=6, res=600)
ggplot(smp.krr, aes(x=percent, y=taxa, colour=tool, shape=min)) + geom_point(size=2) + 
  stat_summary(geom = "line", fun = mean, size=1) + 
  scale_y_continuous(trans='log10') + 
  scale_colour_manual(values=tool_colors, name='Tool',
                      labels=c('Kaiju', 'RAT -mcr', 'RAT -mc')) + 
  scale_shape_discrete(labels=c('All taxa', 'Taxa>=0.001%'), name='') +
  xlab('Percent reads sampled') + ylab('Number of detected taxa') +
  theme(
    panel.background = element_rect(fill='white', colour='#d9d9d9'),
    panel.grid.major = element_line(colour='#bbbbbb'),
    panel.grid.minor = element_line(colour='#dddddd'),
    legend.key.size = unit(4,"mm"),
    text = element_text(size=7)
  )
#dev.off()

ggsave('../../plots/rarefaction/20230428.W23-2.krr.wmin.nomin.svg', width=160, height=80, units='mm')

# insert ggplot code


##### CAMI
data.all<-read.csv('../all.cami.csv')
library(dplyr)

mean_data <- group_by(data.all, sample, tool, percent) %>%
  summarise(taxa = mean(taxa, na.rm = TRUE))

# plot.all<-
  ggplot(mean_data, aes(x=percent, y=taxa, colour=tool)) + 
  geom_line() + geom_point(size=1, alpha=0.2) + facet_wrap(~sample) +     
    theme(panel.background = element_rect(fill='white'),
          panel.grid.major = element_line(colour='#BBBBBB'),
          panel.grid.minor = element_line(colour='#DDDDDD'))
  ggsave('../../plots/rarefaction/cami/all_facet_wrap.png', width=6, height=6)
  




myplots <- vector('list', 10)
i=1

for (s in c('smp6','smp13','smp23','smp25','smp26','smp30','smp33','smp34','smp38','smp53')) {
smp.kaiju<-read.csv(paste(s, 'kaiju.csv', sep='.'))
smp.rat<-read.csv(paste(s, 'rat.csv', sep='.'))
smp.rat.r<-read.csv(paste(s, 'rat_rob.csv', sep='.'))
smp.cami<-read.csv(paste(s, 'cami.csv', sep='.'))
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r, smp.cami)
p.all<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle(s) # + scale_y_continuous(trans='log10')
myplots[[i]]<-p.all
i=i+1
ggsave(paste('../../plots/rarefaction/CAMI/', s, '.all.wmin.png', sep=''), width=5, height=4)

}





smp.kaiju<-read.csv('W19-1.kaiju.nomin.csv')
smp.rat<-read.csv('W19-1.rat.nomin.csv')
smp.rat.r<-read.csv('W19-1.rat_r.nomin.csv')
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r)
p.rat<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle("W19-1") #+ scale_y_continuous(trans='log10')
ggsave('../../plots/rarefaction/W19-1.all.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W19-2.kaiju.nomin.csv')
smp.rat<-read.csv('W19-2.rat.nomin.csv')
smp.rat.r<-read.csv('W19-2.rat_r.nomin.csv')
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r)
p.rat<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle("W19-2") + scale_y_continuous(trans='log10')
ggsave('../plots/rarefaction/W19-2.all.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W19-3.kaiju.nomin.csv')
smp.rat<-read.csv('W19-3.rat.nomin.csv')
smp.rat.r<-read.csv('W19-3.rat_r.nomin.csv')
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r)
p.rat<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle("W19-3") + scale_y_continuous(trans='log10')
ggsave('../plots/rarefaction/W19-3.all.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W19-4.kaiju.nomin.csv')
smp.rat<-read.csv('W19-4.rat.nomin.csv')
smp.rat.r<-read.csv('W19-4.rat_r.nomin.csv')
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r)
p.rat<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle("W19-4") + scale_y_continuous(trans='log10')
ggsave('../plots/rarefaction/W19-4.all.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W19-5.kaiju.nomin.csv')
smp.rat<-read.csv('W19-5.rat.nomin.csv')
smp.rat.r<-read.csv('W19-5.rat_r.nomin.csv')
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r)
p.rat<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle("W19-5") + scale_y_continuous(trans='log10')
ggsave('../plots/rarefaction/W19-5.all.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W19-6.kaiju.nomin.csv')
smp.rat<-read.csv('W19-6.rat.nomin.csv')
smp.rat.r<-read.csv('W19-6.rat_r.nomin.csv')
smp.all<-rbind(smp.kaiju, smp.rat, smp.rat.r)
p.rat<-ggplot(smp.all, aes(x=percent, y=taxa, colour=tool)) + geom_point() + 
  ggtitle("W19-6") + scale_y_continuous(trans='log10')
ggsave('../plots/rarefaction/W19-6.all.nomin.png', width=5, height=4)










smp.kaiju<-read.csv('W19-2.kaiju.nomin.csv')
smp.rat<-read.csv('W19-2.rat.nomin.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W19-2")
ggsave('../plots/rarefaction/W19-2.both.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W19-3.kaiju.nomin.csv')
smp.rat<-read.csv('W19-3.rat.nomin.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W19-3")
ggsave('../plots/rarefaction/W19-3.both.nomin.png', width=5, height=4)


smp.kaiju<-read.csv('W22-1.kaiju.nomin.csv')
smp.rat<-read.csv('W22-1.rat.nomin.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W22-1")
ggsave('../plots/rarefaction/W22-1.both.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W22-3.kaiju.nomin.csv')
smp.rat<-read.csv('W22-3.rat.nomin.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W22-3")
ggsave('../plots/rarefaction/W22-3.both.nomin.png', width=5, height=4)

smp.kaiju<-read.csv('W22-2.kaiju.nomin.csv')
smp.rat<-read.csv('W22-2.rat.nomin.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W22-2")
ggsave('../plots/rarefaction/W22-2.both.nomin.png', width=5, height=4)









smp.kaiju<-read.csv('W19-1.kaiju.min0.001.csv')
smp.rat<-read.csv('W19-1.rat.min0.001.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W19-1")
ggsave('../plots/rarefaction/W19-1.both.min0.001.png', width=5, height=4)

smp.kaiju<-read.csv('W19-2.kaiju.min0.001.csv')
smp.rat<-read.csv('W19-2.rat.min0.001.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W19-2")
ggsave('../plots/rarefaction/W19-2.both.min0.001.png', width=5, height=4)

smp.kaiju<-read.csv('W19-3.kaiju.min0.001.csv')
smp.rat<-read.csv('W19-3.rat.min0.001.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W19-3")
ggsave('../plots/rarefaction/W19-3.both.min0.001.png', width=5, height=4)


smp.kaiju<-read.csv('W22-1.kaiju.min0.001.csv')
smp.rat<-read.csv('W22-1.rat.min0.001.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W22-1")
ggsave('../plots/rarefaction/W22-1.both.min0.001.png', width=5, height=4)

smp.kaiju<-read.csv('W22-3.kaiju.min0.001.csv')
smp.rat<-read.csv('W22-3.rat.min0.001.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W22-3") +
  scale_y_continuous(trans="log10")
ggsave('../plots/rarefaction/W22-3.both.min0.001.png', width=5, height=4)

smp.kaiju<-read.csv('W22-2.kaiju.min0.001.csv')
smp.rat<-read.csv('W22-2.rat.min0.001.csv')
smp.both<-rbind(smp.kaiju, smp.rat)
p.rat<-ggplot(smp.both, aes(x=percent, y=taxa, colour=tool)) + geom_point() + ggtitle("W22-2")
ggsave('../plots/rarefaction/W22-2.both.min0.001.png', width=5, height=4)

