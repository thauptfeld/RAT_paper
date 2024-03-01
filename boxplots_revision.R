
setwd("C:/Users/acer/OneDrive - Flinders/RATrevision/scripts")
setwd("C:/Users/haup0007/OneDrive - Flinders/RATrevision/scripts")
library(ggplot2)
library(reshape2)
library(Hmisc)
library(dplyr)
library(ape)
library(vegan)


#Algorithm Stats - NCBI
tool_colors=c('#309f74','#b9f291', '#2278b4', '#7fb1d3',
              '#f6e5ad', "#f72349", '#5e5959')
stats_ncbi <- read.csv('../evaluations/20231020.CAMI2_marine.csv', header=F)

stats_ncbi$V4 <- factor(stats_ncbi$V4, levels=c("diamond","kaiju", "centrifuge", "kraken",
                                                "r_sensitive",'r_contig', 'r_robust'))
tool_labels=c("DIAMOND","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr (CAMI2)", "RAT -mcr (MetaBAT2)", "RAT -mc")

# tool_colors=c('#309f74',
#               '#f6dc85', "#f72349", '#5e5959')
# stats_ncbi <- read.csv('../revision/evaluations/20231009.CAMI2_marine_onlyRAT.csv', header=F)
# stats_ncbi$V4 <- factor(stats_ncbi$V4, levels=c("diamond",
#                                                       "r_sensitive", 'r_contig', 'r_robust'))
# tool_labels=c("Diamond",
#               "RAT CAMI genomes", "RAT MetaBAT2 MAGs", "RAT without MAGs", "RAT robust")




p <- ggplot(stats_ncbi, aes(x=V2, y=V3, fill=V4)) +
  geom_violin(scale = "width", width=1, color=NA) +facet_wrap(~V1)+ ylab("Fraction of Reads") + xlab("Taxonomic Rank") + 
  scale_fill_manual(name = "Tool", values=tool_colors,
                    labels = tool_labels) +  
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", 
  #              position=position_dodge(1), size=0.1) +
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
ggsave('../plots/20231129_marine_correct_incorrect_onlyRAT.png', height=4, width=7.7, dpi=600)


#Algorithm Stats - Plant
tool_colors=c('#309f74','#b9f291', '#2278b4', '#7fb1d3',
              '#f6e5ad', "#f72349", '#5e5959')
stats_ncbi <- read.csv('../evaluations/20231020_CAMI2_plant.csv', header=F)

stats_ncbi$V4 <- factor(stats_ncbi$V4, levels=c("diamond","kaiju", "centrifuge", "kraken",
                                                "r_sensitive",'r_contig', 'r_robust'))
tool_labels=c("Diamond","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr (CAMI2)", "RAT -mcr (MetaBAT2)", "RAT -mc")

# tool_colors=c('#309f74',
#               '#f6dc85', "#f72349", '#5e5959')
# stats_ncbi <- read.csv('../revision/evaluations/20231009.CAMI2_marine_onlyRAT.csv', header=F)
# stats_ncbi$V4 <- factor(stats_ncbi$V4, levels=c("diamond",
#                                                       "r_sensitive", 'r_contig', 'r_robust'))
# tool_labels=c("Diamond",
#               "RAT CAMI genomes", "RAT MetaBAT2 MAGs", "RAT without MAGs", "RAT robust")




p <- ggplot(stats_ncbi, aes(x=V2, y=V3, fill=V4)) +
  geom_violin(scale = "width", width=1, color=NA) +facet_wrap(~V1)+ ylab("Fraction of Reads") + xlab("Taxonomic Rank") + 
  scale_fill_manual(name = "Tool", values=tool_colors,
                    labels = tool_labels) +  
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", 
  #              position=position_dodge(1), size=0.01) +
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
ggsave('../plots/20231031_plant_correct_incorrect_onlyRAT.png', height=4, width=7.7, dpi=600)








# Read annotation barplot
# # stats <- read.csv('../benchmark/mousegut_fraction.csv')
# stats <- read.csv('../RAT_wbin_reads_per_step_frac.csv')
# 
# stats.m <- melt(stats)
# stats.m$X <- factor(stats.m$X, levels=c('unclassified', "read_dm", "contig_dm", "contig", 'bin'))
# p <- ggplot(stats.m, aes(x=variable, y=value, fill=X)) + geom_bar(stat='identity') +
#   scale_fill_manual(values=c('#7fb1d3','#f6dc85', "#f72349", '#ff845e', '#5e5959'), name='Step') +
#   ylab("Fraction of reads classified") +
#   xlab("sample") +
#   theme(
# axis.text.x = element_text(angle=270, hjust=0.0, vjust=0.4), 
# panel.background = element_rect(fill='#FFFFFF'), 
# panel.grid.major = element_line(colour='#BBBBBB'), 
# panel.grid.minor =  element_line(colour='#CCCCCC')
# )
# 
# 
# 
# p
# ggsave('20230227_cami_reads_per_step.png', height=4, width=6, dpi=600)



# Tool_stats
elfcols=c('#364B9A', '#990099', '#DD00AA', '#DD00AA', '#98CAE1',
          '#98CAE1', '#98CAE1', '#98CAE1', '#98CAE1', '#98CAE1',
          '#98CAE1', '#009900', '#EAECCC', '#FEDA8B', '#FDB366',
          '#555555', '#DD3D2D', '#DD3D2D', '#A50026')

stats <- read.csv('../evaluations/20231024.all_tools.plant.0.001.csv', header=T)

tool_colors=c('#50bf94','#b9f291', '#2278b4', '#7fb1d3', 
              '#f6e5ad', "#f72349", '#5e5959')

stats$tool <- factor(stats$tool, levels=c("bracken","kaiju", "centrifuge", "kraken2",
                                                "sensitive", 'contig', 'robust'))
tool_labels=c("Bracken","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr", "RAT -cr", "RAT -mc")

#stats2=stats[which(stats$tool=='RAT' | stats$tool=='CK_v0'),]
stats2=stats[which(stats$metric=='l1'),]
stats3=stats2[,-c(3)]


ggplot(stats3, aes(x=rank, y=value, fill=tool)) +
  geom_violin(scale = "width", width=1, color=NA) + ylab("L1 Distance") + xlab("Taxonomic Rank")  + # facet_wrap(~metric, scales='free') +
  scale_x_discrete(limits=c("phylum", "class", "order", "family", "genus", "species")) +
  scale_fill_manual(values=tool_colors, name="Tool", labels=tool_labels) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", 
               position=position_dodge(1), size=0.1) +
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=270, hjust=0.0, vjust=0.4),
    panel.border = element_rect(color='#d9d9d9', fill=NA, size=1),
    panel.background = element_rect(fill='#FFFFFF'), 
    legend.key.size = unit(2,"line"),
    panel.grid.major = element_blank(), #element_line(colour='#BBBBBB'), 
    panel.grid.minor =  element_blank() #element_line(colour='#CCCCCC')
  )


ggsave('../plots/20231031_L1_CAMI2_plant_all_tools.png', height=7, width=8, dpi=600)



### Unifrac
unifrac=read.csv('../evaluations/20231028.marine_unifrac.csv')
ggplot(unifrac, aes(x=tool, y=value, fill=tool)) + geom_boxplot() + ylab('Unifrac distance') + xlab("Tool") +
  scale_fill_manual(values=elfcols) + scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10)) +
  theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'), panel.grid.minor =  element_line(colour='#CCCCCC')
  )
ggsave('Unifrac_CAMI1_all_tools.png', height=6, width=6.5, dpi=600)



sens.prec.=read.csv('../evaluations/20231024.CAMI2_plant.sens.prec.csv')
sens.prec=sens.prec.[which(sens.prec.$rank!='superkingdom'),]

tool_colors=c('#50bf94','#b9f291', '#2278b4', '#7fb1d3', 
              '#f6e5ad', "#f72349", '#5e5959')

sens.prec$tool <- factor(sens.prec$tool, levels=c("bracken","kaiju", "centrifuge", "kraken2",
                                          "sensitive", 'contig', 'robust'))
tool_labels=c("Bracken","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr", "RAT -cr", "RAT -mc")

tiff("../plots/20231031_sens_prec_CAMI2_plant_all_tools_0.001_BIG.tiff", units="in", width=9.5, height=7, res=600)

sens.prec$rank_f = factor(sens.prec$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(sens.prec, aes(x=sensitivity, y=precision, color=tool, shape=rank_f)) + 
  geom_point(alpha=0.8, size=4) + xlab("Sensitivity") + ylab("Precision") +
  scale_color_manual(values=tool_colors, labels=tool_labels, name="Tool") +
  scale_shape_manual(values=c(3, 6, 15,16,17,18), name="Rank") +
  scale_x_continuous(limits=c(0,1)) +
  theme(
    panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'),
    panel.grid.minor =  element_line(colour='#CCCCCC'),
    text=element_text(size=20), 
    legend.key.size = unit(1.8, 'line')
  )

dev.off()




sens.prec.=read.csv('../evaluations/20231024.CAMI2_marine.sens.prec.csv')
sens.prec=sens.prec.[which(sens.prec.$rank!='superkingdom'),]

tool_colors=c('#50bf94','#b9f291', '#2278b4', '#7fb1d3', 
              '#f6e5ad', "#f72349", '#5e5959')

sens.prec$tool <- factor(sens.prec$tool, levels=c("bracken","kaiju", "centrifuge", "kraken2",
                                                  "sensitive", 'contig', 'robust'))
tool_labels=c("Bracken","Kaiju", "Centrifuge", "Kraken2",
              "RAT -mcr", "RAT -cr", "RAT -mc")
sens.prec$rank_f = factor(sens.prec$rank, levels=c("phylum", "class", "order", "family", "genus", "species"))
ggplot(sens.prec, aes(x=sensitivity, y=precision, color=tool, shape=rank_f)) + 
  geom_point(alpha=0.8, size=4) + xlab("Sensitivity") + ylab("Precision") +
  ylim(0,1) +
  scale_color_manual(values=tool_colors, labels=tool_labels, name="Tool") +
  scale_shape_manual(values=c(3, 6, 15,16,17,18), name="Rank") +
  scale_x_continuous(limits=c(0,1)) +
  theme(
    panel.background = element_rect(fill='#FFFFFF'), 
    panel.grid.major = element_line(colour='#BBBBBB'),
    panel.grid.minor =  element_line(colour='#CCCCCC'),
    text=element_text(size=20), 
    legend.key.size = unit(1.8, 'line')
  )
ggsave('../plots/20231024_sens_prec_CAMI2_marine_all_tools_0.001_small.png', height=7, width=8, dpi=600)




### Memory stuff



### Bar Chart taxa

