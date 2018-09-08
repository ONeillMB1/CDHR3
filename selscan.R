require(ggplot2)
library(viridis)

sel <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/170405_selscan_totpop_maf05_pel_and_chb_maf04.txt", header=T, sep = "\t")


sel$pop <- factor(sel$pop, levels=rev(unique((sel$pop)[order(sel$ord)])))
#ord <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/1k_pop.txt", header=F, sep="\t")
#names(ord) <- c("pop", "pop_desc", "superpop", "ord")

#levels(sel$pop) <- rev(levels(factor(ord$pop, levels=unique(ord$pop))))
#levels(sel$sex) <- c("None", "female", "male")

ggplot(sel) + geom_boxplot(aes(x=stat, y = per_rank))

#ggplot(sel) + geom_boxplot(aes(x=sex, y = per_rank)) + facet_wrap(~stat)



ssp <- ggplot(sel, aes(y=pop, x=stat)) + 
  geom_tile(aes(fill=per_rank)) + 
  geom_text(aes(label = paste(signif(val, 2), " (", signif(per_rank, 2), "%)", sep="")), size=3) +
  #facet_wrap(~sex) +
  scale_fill_gradientn(colors=c("#FFFF66","#FF3300"), 
                       limits=c(0,100), 
                       values=c(0.0, .90, .95, 1.0), 
                       breaks=c(0,50, 90, 95, 100)) +
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="top",
        legend.title = element_blank(),
        legend.key.width = unit(0.41, "in")) +
  xlab("") +
  ylab("")
ssp  

ggsave(plot=ssp, filename = "C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/170405_totpop.pdf", width=4, height=9, units="in", dpi=300 )
