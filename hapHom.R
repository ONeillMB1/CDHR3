require(ggplot2)

dat <- read.table("C:/Users/Mary/PepLab/data/CDHR3/chr7_105000000-106000000_maf0.1.txt", header=F, sep="\t", na.strings="NA")

dat.p <- dat[!is.na(dat$V3),]


dat.p$length <- dat.p$V4 - dat.p$V3
dat.p <- dat.p[order(dat.p$V2, dat.p$length),]
dat.p$order <- seq(1,length(dat.p$V1))

pop <- read.table("C:/Users/Mary/PepLab/data/CDHR3/integrated_call_samples_v3.20130502.ALL.panel")


dat.p$pop <- pop[match(dat.p$V1, pop$V1), 'V2']
dat.p$superpop <- pop[match(dat.p$V1, pop$V1), 'V3']
dat.p$sex <- pop[match(dat.p$V1, pop$V1), 'V4']




dat.ord <- dat.p[order(dat.p$V2, dat.p$superpop, dat.p$pop, dat.p$length),]
dat.ord$order <- seq(1, length(dat.ord$V1))


ggplot(dat.ord) + 
  #facet_wrap(~sex, scales = "free_y") +
  geom_segment(aes(x=V3, y=order, xend=V4, yend=order, color=V2)) +
  xlab("Chromosome 7") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
  #scale_x_continuous(limits=c(0e06, 1e06)) +

GG <- dat.p[as.character(dat.p$V2) == "GG", 'length']
AA <- dat.p[as.character(dat.p$V2) == "AA", 'length']





aggregate(length ~ sex + V2 + pop, data = dat.p, FUN=mean)
