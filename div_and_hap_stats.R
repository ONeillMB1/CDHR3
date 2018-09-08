require(ggplot2)

GWD <- read.table("C:/Users/Mary/PepLab/data/CDHR3/GWD_test_iHS.norm", header=F, sep = "\t")

GWD_0 <- GWD[(GWD$V8 == 0 & GWD$V9 == 0),]
GWD_1 <- GWD[(GWD$V8 != 0 & GWD$V9 != 0),]

ggplot(GWD_1) + geom_point(aes(x=V2, y=V7))
ggplot(GWD_0) + geom_point(aes(x=V2, y=abs(V7)))

GWD[GWD$V7 < -4,]



####

ehh <- read.table("C:/Users/Mary/PepLab/data/CDHR3/GWD_test_ehh", header=F, sep = "\t")

ggplot(ehh) + geom_line(aes(x=V1, y=V2, colour=allele))


ehh2 <- read.table("C:/Users/Mary/PepLab/data/CDHR3/test_ehh", header=F, sep = "\t")
ggplot(ehh2) + geom_line(aes(x=V1, y=V2, colour=as.factor(V3)))



####

all <- read.table("C:/Users/Mary/PepLab/data/CDHR3/thousand_genomes_meltEHH.txt", header=F, sep = "\t")

names(all) <- c("pop", "pos", "EHH", "allele", "relpos")
all$allele <- as.factor(all$allel)
all$relpos <- as.factor(all$relpos)
str(all)


ggplot(all) + 
  facet_wrap(~pop) +
  geom_line(aes(pos, EHH, color=allele)) +
  theme_bw() +
  xlab("Genomic Position (MB)") +
  ylab("Extended Haplotype Heterozygosity (EHH)") +
  scale_x_continuous(limits=c(105500000,105800000), breaks=c(105500000,105600000,105700000,105800000), labels=c("105.5", "105.6", "105.7", "105.8")) +
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        legend.title=element_blank(),
        axis.title.x=element_text(),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(hjust = 0.5, vjust = 1, size=rel(1)),
        #axis.title.y=element_text(hjust = 0.5, vjust = 0.05),
        axis.line=element_line()#,
        #strip.text.x=element_blank(),
        #strip.background=element_blank(),
        #plot.margin = unit(c(0.05, 0.5, 0.05, 0.5), "lines")
  )


###
all2 <- read.table("C:/Users/Mary/PepLab/data/CDHR3/thousand_genomes_meltEHH_0.05.txt", header=F, sep = "\t")

names(all2) <- c("pop", "pos", "EHH", "allele", "relpos")
all2$allele <- as.factor(all2$allel)
all2$relpos <- as.factor(all2$relpos)
str(all2)


ggplot(all2) + 
  facet_wrap(~pop) +
  geom_line(aes(pos, EHH, color=allele)) +
  xlab("Genomic Position") +
  ylab("Extended Haplotype Heterozygosity (EHH)") +
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        legend.title=element_blank(),
        axis.title.x=element_text(),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(hjust = 0.5, vjust = 1, size=rel(1)),
        #axis.title.y=element_text(hjust = 0.5, vjust = 0.05),
        axis.line=element_line()#,
        #strip.text.x=element_blank(),
        #strip.background=element_blank(),
        #plot.margin = unit(c(0.05, 0.5, 0.05, 0.5), "lines")
  )

###
TD <- read.table("C:/Users/Mary/PepLab/data/CDHR3/thousand_genomes_TajimaD_10K.txt", header=F, sep = "\t")
names(TD) <- c("POP", "CHROM", "BIN_START", "N_SNPS", "TajimaD")
str(TD)

td.p <- ggplot(TD) + 
  facet_wrap(~POP) +
  geom_line(aes((BIN_START+25000), TajimaD)) +
  geom_vline(xintercept=105658451, colour="red") +
  xlab("Genomic Position") +
  ylab("Tajima's D") +
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        legend.title=element_blank(),
        axis.title.x=element_text(),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(hjust = 0.5, vjust = 1, size=rel(1)),
        #axis.title.y=element_text(hjust = 0.5, vjust = 0.05),
        axis.line=element_line()#,
        #strip.text.x=element_blank(),
        #strip.background=element_blank(),
        #plot.margin = unit(c(0.05, 0.5, 0.05, 0.5), "lines")
  )

td.p + scale_x_continuous(limits=c(105203657,106076877))

td.p + scale_x_continuous(limits=c(105403657,105876877))

td.p + scale_x_continuous(limits=c(105503657,105776877))

td.reg <- TD[TD$BIN_START < 105658541 & (TD$BIN_START+50000) > 105658451,]

##
PI <- read.table("C:/Users/Mary/PepLab/data/CDHR3/thousand_genomes_Pi_w10K_s1K.txt", header=F, sep = "\t")
names(PI) <- c("POP", "CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI")
str(PI)
PI$x <- (PI$BIN_END+PI$BIN_START)/2

pi.p <- ggplot(PI) + 
  facet_wrap(~POP) +
  geom_line(aes(x, PI)) +
  geom_vline(xintercept=105658451, colour="red") +
  xlab("Genomic Position") +
  ylab(expression(pi)) +
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        legend.title=element_blank(),
        axis.title.x=element_text(),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(hjust = 0.5, vjust = 1, size=rel(1)),
        #axis.title.y=element_text(hjust = 0.5, vjust = 0.05),
        axis.line=element_line()#,
        #strip.text.x=element_blank(),
        #strip.background=element_blank(),
        #plot.margin = unit(c(0.05, 0.5, 0.05, 0.5), "lines")
  )

pi.p + scale_x_continuous(limits=c(105203657,106076877))

pi.p + scale_x_continuous(limits=c(105403657,105876877))

pi.p + scale_x_continuous(limits=c(105503657,105776877))

### Nuc Diversity & EHH (vcflib, default 20 SNPs)

div <- read.table("C:/Users/Mary/PepLab/data/CDHR3/thousand_genomes_seqDiv.txt", header=F, sep = "\t")

names(div) <- c("pop", "chrom", "start", "end", "pi", "ehh")
str(div)
div$x <- (div$end+div$start)/2

ggplot(div) + 
  facet_wrap(~pop) +
  geom_line(aes(x, pi)) +
  geom_vline(xintercept=105658451, colour="red")


ggplot(div) + 
  facet_wrap(~pop) + 
  geom_point(aes(x, ehh)) +
  geom_vline(xintercept=105658451, colour="red")


#iHS

ihs <- read.table("C:/Users/Mary/PepLab/data/CDHR3/thousand_genomes_iHS_0.05_norm.txt", header=F, sep = "\t")
names(ihs) <- c("pop", "chrom", "pos", "alleleFreq", "EHH_alt", "EHH_ref", "iHS", "norm_iHS", "intFail_1", "intFail_2")

ihs_0 <- ihs[(ihs$intFail_1 == 0 & ihs$intFail_2 == 0),]
ihs_1 <- ihs[(ihs$intFail_1 != 0 & ihs$intFail_2 != 0),]

neg <- ihs_0[ihs_0$iHS < 0,]

ggplot(neg) + 
  facet_wrap(~pop) + 
  geom_point(aes(pos, iHS), size = .5) +
  geom_vline(xintercept=105658451, colour="red")

ggplot(ihs_0) + 
  facet_wrap(~pop) + 
  geom_point(aes(pos, abs(iHS)), size = .05) +
  geom_vline(xintercept=105658451, colour="red") +
  xlab("Genomic Position") +
  ylab("Integrated Haplotype Score (iHS)") +
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        legend.title=element_blank(),
        axis.title.x=element_text(),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(hjust = 0.5, vjust = 1, size=rel(1)),
        #axis.title.y=element_text(hjust = 0.5, vjust = 0.05),
        axis.line=element_line()#,
        #strip.text.x=element_blank(),
        #strip.background=element_blank(),
        #plot.margin = unit(c(0.05, 0.5, 0.05, 0.5), "lines")
  )








#------
  
age <- function(p) {
  E = ((-2*p)/(1-p)) * log(p)
  return(E)
}


age(0.02)
