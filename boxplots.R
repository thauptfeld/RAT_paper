
setwd("~/PhD/RAT/plots")
library(ggplot2)
library(reshape2)
library(Hmisc)
library(dplyr)
library(ape)
library(vegan)
library(svglite)


#Algorithm Stats - NCBI
# tool_colors=c('#309f74','#b9f291', '#2278b4', '#7fb1d3',
#               '#f6dc85', '#ff845e', "#f72349", '#5e5959')
# stats_ncbi <- read.csv('../stats_CAMI2/20230127_CAMI2_w_unclassified_CORRECT.csv')
# 
# stats_ncbi$group <- factor(stats_ncbi$group, levels=c("diamond","kaiju", "centrifuge", "kraken",
#                                                 "r_wbin", "r_wmetabat", 'r_nobin', 'r_robust'))
# tool_labels=c("Diamond","Kaiju", "Centrifuge", "Kraken2",
#               "RAT_CAMI", "RAT_metabat", "RAT_nobin", "RAT_robust")

tool_colors=c('#309f74',
              '#f6dc85', '#ff845e', "#f72349", '#5e5959')
stats_ncbi <- read.csv('../stats_CAMI2/20220902_CAMI2_w_unclassified_CORRECT_onlyRAT.csv')
stats_ncbi$group <- factor(stats_ncbi$group, levels=c("diamond",
                                                      "r_wbin", "r_wmetabat", 'r_nobin', 'r_robust'))
tool_labels=c("Diamond",
              "RAT CAMI genomes", "RAT MetaBAT2 MAGs", "RAT without MAGs", "RAT robust")




p <- ggplot(stats_ncbi, aes(x=rank, y=value, fill=group)) +
  geom_violin(scale = "width", width=1, color=NA) +facet_wrap(~measure)+ ylab("Fraction of Reads") + xlab("Taxonomic Rank") + 
  scale_fill_manual(name = "Tool", values=tool_colors,
                    labels = tool_labels) +  
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", 
               position=position_dodge(1), size=0.1) +
  scale_x_discrete(limits=c("phylum", "class", "order", "family", "genus", "species")) +
  scale_y_continuous(limits=c(-0.1, 1)) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0.0, vjust=0.4), 
    panel.background = element_rect(fill='#FFFFFF', colour="#d9d9d9"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor =  element_blank(),
    text = element_text(size=10)
  )

p
ggsave('../plots/20230222_correct_incorrect_wuncl_wdiamond_CORRECT_onlyRAT.png', height=4, width=7.7, dpi=600)






tool_colors=c('#309f74','#b9f291', '#2278b4', '#7fb1d3', 
              '#f6dc85', '#ff845e', "#f72349", '#5e5959')
stats_ncbi <- read.csv('../stats_CAMI2/20220414_CAMI2_TPR.csv')
stats_ncbi$group <- factor(stats_ncbi$group, levels=c("diamond","kaiju", "centrifuge", "kraken",
                                                      "r_wbin", "r_wmetabat", 'r_nobin', 'r_robust'))
tool_labels=c("Diamond","Kaiju", "Centrifuge", "Kraken2",
              "RAT_genomebins", "RAT_metabat", "RAT_nobin", "RAT_robust")



p <- ggplot(stats_ncbi, aes(x=rank, y=TPR, fill=group)) +
  geom_violin(scale = "width", width=1) + ylab("True Positive Rate") + xlab("Taxonomic Rank") + 
  scale_fill_manual(name = "Tool", values=tool_colors,
                    labels = tool_labels) +  
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", 
               position=position_dodge(1), size=0.1) +
  scale_x_discrete(limits=c("phylum", "class", "order", "family", "genus", "species")) +
  scale_y_continuous(limits=c(-0.1, 1)) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0.0, vjust=0.4), 
    panel.background = element_rect(fill='#FFFFFF', colour='#d9d9d9'), 
    panel.grid.major = element_blank(), 
    panel.grid.minor =  element_blank(),
    text = element_text(size=10)
  )

p
ggsave('../plots/20220414_TPR.png', height=4, width=4.5, dpi=600)


p <-ggplot(stats_ncbi, aes(x=unclassified, y=TPR, colour=group, shape=rank)) + geom_point() +
  scale_colour_manual(values=tool_colors)
p



# Read annotation barplot
# stats <- read.csv('../benchmark/mousegut_fraction.csv')
stats <- read.csv('../RAT_wbin_reads_per_step_frac.csv')

stats.m <- melt(stats)
stats.m$X <- factor(stats.m$X, levels=c('unclassified', "read_dm", "contig_dm", "contig", 'bin'))
p <- ggplot(stats.m, aes(x=variable, y=value, fill=X)) + geom_bar(stat='identity') +
  scale_fill_manual(values=c('#7fb1d3','#f6dc85', "#f72349", '#ff845e', '#5e5959'), name='Step') +
  ylab("Fraction of reads classified") +
  xlab("sample") +
  theme(
axis.text.x = element_text(angle=270, hjust=0.0, vjust=0.4), 
panel.background = element_rect(fill='#FFFFFF'), 
panel.grid.major = element_line(colour='#BBBBBB'), 
panel.grid.minor =  element_line(colour='#CCCCCC')
)



p
ggsave('20230227_cami_reads_per_step.png', height=4, width=6, dpi=600)



# Tool_stats
elfcols=c('#364B9A', '#990099', '#DD00AA', '#DD00AA', '#98CAE1', 
          '#98CAE1', '#98CAE1', '#98CAE1', '#98CAE1', '#98CAE1',
          '#98CAE1', '#009900', '#EAECCC', '#FEDA8B', '#FDB366',
          '#555555', '#DD3D2D', '#DD3D2D', '#A50026')

stats <- read.csv('../stats_CAMI2/20230321_all_tools_stats_0.001.csv')

tool_colors=c('#50bf94','#b9f291', '#2278b4', '#7fb1d3', 
              '#f6dc85', '#ff845e', "#f72349", '#5e5959')

stats$tool <- factor(stats$tool, levels=c("bracken","kaiju", "centrifuge", "kraken2",
                                                "wbin", "wmetabat", 'nobin', 'robust'))
tool_labels=c("Bracken","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr (CAMI)", "RAT -mcr (MetaBAT2)", "RAT -cr", "RAT -mc")

#stats2=stats[which(stats$tool=='RAT' | stats$tool=='CK_v0'),]
stats2=stats[which(stats$metric=='l1'),]
stats3=stats2[,-c(3)]

# tiff("20230503_L1_CAMI2_all_tools_BIG.tiff", units="in", width=9, height=7, res=600)

ggplot(stats3, aes(x=rank, y=value, fill=tool)) +
  geom_violin(scale = "width", width=1, color=NA) + ylab("L1 Distance") + xlab("Taxonomic Rank")  + # facet_wrap(~metric, scales='free') +
  scale_x_discrete(limits=c("phylum", "class", "order", "family", "genus", "species")) +
  scale_fill_manual(values=tool_colors, name="Tool", labels=tool_labels) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", 
               position=position_dodge(1), size=0.05) +
  theme(
    text = element_text(size=9),
    axis.text.x = element_text(angle=270, hjust=0.0, vjust=0.4),
    panel.border = element_rect(color='#d9d9d9', fill=NA, size=1),
    panel.background = element_rect(fill='#FFFFFF'), 
    legend.key.size = unit(4,"mm"),
    legend.box.margin = margin(0,-5,0,-10),
    panel.grid.major = element_blank(), #element_line(colour='#BBBBBB'), 
    panel.grid.minor =  element_blank() #element_line(colour='#CCCCCC')
  )

# dev.off()

ggsave('../plots/20240228_L1_CAMI2_all_tools.svg', height=90, width=130, units = 'mm')

unifrac=stats[which(stats$metric=='unifrac'),]
ggplot(unifrac, aes(x=tool, y=value, fill=tool)) + geom_boxplot() + ylab('Unifrac distance') + xlab("Tool") +
  scale_fill_manual(values=elfcols) + scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10)) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'), panel.grid.minor =  element_line(colour='#CCCCCC')
  )
ggsave('Unifrac_CAMI1_all_tools.png', height=6, width=6.5, dpi=600)


# L1 Facet Wrap
l1=stats[which(stats$metric=='l1'),]
l1$rank_f = factor(l1$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(l1, aes(x=tool, y=value, fill=tool)) + geom_boxplot() + ylab('L1 Norm') + xlab("Tool") +
  scale_fill_manual(values=elfcols) +
  facet_wrap(~rank_f) + theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'), panel.grid.minor =  element_line(colour='#CCCCCC')
  )
  
ggsave('l1_CAMI1_all_tools.png', height=12, width=13, dpi=600)


# L1 Line plot
l1=read.csv("l1_avg.csv")
l1$rank_f = 10
l1$rank_f[l1$rank=="class"]<-20
l1$rank_f[l1$rank=="order"]<-30
l1$rank_f[l1$rank=="family"]<-40
l1$rank_f[l1$rank=="genus"]<-50
l1$rank_f[l1$rank=="species"]<-60
l1$thickness='thin'
l1$thickness[l1$tool=="RAT"]<-'thick'
ggplot(l1, aes(x=rank_f, y=avg, color=tool, group=tool, size=thickness)) + geom_line() + ylab('L1 Norm') + xlab("Tool") +
  scale_x_continuous(limits=c(10,60), breaks=c(10,20,30,40,50,60), labels=c("phylum", "class", "order", "family", "genus", "species"))+
  scale_color_manual(values=elfcols) + geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), size = 0.3, width=12, color="black", position=position_dodge(0.8)) + 
  scale_size_manual(values=c(1.5, 0.5)) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'), panel.grid.minor =  element_line(colour='#CCCCCC')
  )

ggsave('l1_CAMI1_lineplot.png', height=6, width=6, dpi=600)



sens.prec.=read.csv('../stats_CAMI2/20230321_all_tools_stats_0.001_sens_prec.csv')
sens.prec=sens.prec.[which(sens.prec.$rank!='superkingdom'),]

tool_colors=c('#50bf94','#b9f291', '#2278b4', '#7fb1d3', 
              '#f6dc85', '#ff845e', "#f72349", '#5e5959')

sens.prec$tool <- factor(sens.prec$tool, levels=c("bracken","kaiju", "centrifuge", "kraken2",
                                          "wbin", "wmetabat", 'nobin', 'robust'))
tool_labels=c("Bracken","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr (CAMI)", "RAT -mcr (MetaBAT2)", "RAT -cr", "RAT -mc")

# tiff("20230503_sens_prec_CAMI2_all_tools_0.001_BIG.tiff", units="in", width=9.5, height=7, res=600)

sens.prec$rank_f = factor(sens.prec$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(sens.prec, aes(x=sensitivity, y=precision, color=tool, shape=rank_f)) + 
  geom_point(alpha=0.8, size=2) + xlab("Sensitivity") + ylab("Precision") +
  scale_color_manual(values=tool_colors, labels=tool_labels, name="Tool") +
  scale_shape_manual(values=c(3, 6, 15,16,17,18), name="Rank") +
  scale_x_continuous(limits=c(0,1)) +
  theme(
    panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'),
    panel.grid.minor =  element_line(colour='#CCCCCC'),
    text=element_text(size=7), 
    legend.key.size = unit(3, 'mm'),
    legend.box.margin = margin(0,-5,0,-10)
  )

#dev.off()

ggsave('20240228_sens_prec_CAMI2_all_tools_0.001.svg', height=75, width=88, units='mm')



TP.FP=read.csv('all_tools_CAMI2_TP_FP_0.1.csv')

TP.FP$rank_f = factor(TP.FP$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(TP.FP, aes(x=FP, y=TP, color=tool, shape=rank_f)) + geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=tool_colors) +
  scale_shape_manual(values=c(3, 6, 15,16,17,18)) +
  scale_x_continuous(trans='log10') +
  theme(
    panel.background = element_rect(fill='#FFFFFF'), panel.grid.major = element_line(colour='#BBBBBB'),
    panel.grid.minor =  element_line(colour='#CCCCCC') 
  )

ggsave('TP_FP_CAMI1_all_tools_0.1.png', height=7.5, width=8, dpi=600)





precision=stats[which(stats$metric=='precision'),]
precision$rank_f = factor(precision$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(precision, aes(x=tool, y=value, fill=tool)) + geom_boxplot() + ylab('Precision') + xlab("Tool") +
  facet_wrap(~rank_f) + # scale_x_discrete(limits=c("phylum", "class", "order", "family", "genus", "species")) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'), panel.grid.minor =  element_line(colour='#CCCCCC')
  )

ggsave('precision_CAMI1_all_tools.png', height=12, width=13, dpi=600)

sensitivity=stats[which(stats$metric=='sensitivity'),]
sensitivity$rank_f = factor(sensitivity$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(sensitivity, aes(x=tool, y=value, fill=tool)) + geom_boxplot() + ylab('Sensitivity') + xlab("Tool") +
  facet_wrap(~rank_f) + # scale_x_discrete(limits=c("phylum", "class", "order", "family", "genus", "species")) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'), panel.grid.minor =  element_line(colour='#CCCCCC')
  )

ggsave('sensitivity_CAMI1_all_tools.png', height=12, width=13, dpi=600)






### Memory stuff



### Bar Chart taxa

