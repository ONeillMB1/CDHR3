
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(cowplot)

nsl <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/samp_master_normalized.txt", header=T, sep = "\t", stringsAsFactors = TRUE)
class(nsl)

nsl <- tbl_df(nsl)
class(nsl)

glimpse(nsl)
head(nsl)

nsl <- nsl %>%
  group_by(pop) %>%
  mutate(percrank=rank(abs(norm_nsl))/length(norm_nsl))


aa <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/samp_AA_info.txt", header=T)

aa$derived <- sapply(strsplit(as.character(aa$AA), "|"), `[`, 1)
aa$derived <- as.factor(aa$derived)

#creat unique identifier
aa$snpID <- paste(aa$CHROM, aa$POS, sep="_")
nsl$snpID <- paste(str_replace(as.character(nsl$chr), "chr", ""), nsl$phys_pos, sep="_")




aa$der_AF <- ifelse(as.character(aa$derived) == as.character(aa$ALT), 1-aa$AF, aa$AF)


nslAA <- left_join(nsl, aa, by="snpID")
  
#extract SNPs that have an A,T,C,G in the ancestral allele  
der.m <- nslAA[nslAA$derived == "A" | nslAA$derived == "T" | nslAA$derived == "C" | nslAA$derived == "G",]

der.m$pop_der_AF <- ifelse(as.character(der.m$derived) == as.character(der.m$ALT), 1-der.m$X1_freq, der.m$X1_freq)

pop26 <- levels(der.m$pop)
pop26 <- pop26[c(1, 8:14, 18, 22:35, 39:41)]

pops <- filter(der.m, pop %in% pop26)


meanDAF <- pops %>%
  group_by(locus) %>%
  summarise(meanDAF = mean(pop_der_AF), sdDAF = sd(pop_der_AF))




#transform data to wide
per <- dcast(der.m, locus + snpID + REF + ALT + derived + AF + der_AF + EAS_AF + AMR_AF + AFR_AF + EUR_AF + SAS_AF ~ pop, value.var = c("percrank"))

val <- dcast(der.m, locus + snpID + REF + ALT + derived + AF + der_AF + EAS_AF + AMR_AF + AFR_AF + EUR_AF + SAS_AF ~ pop, value.var = c("norm_nsl"))

freq <- dcast(der.m, locus + snpID + REF + ALT + derived + AF + der_AF + EAS_AF + AMR_AF + AFR_AF + EUR_AF + SAS_AF ~ pop, value.var = c("X1_freq"))

#Summarize

cdhr3 <- nsl[nsl$locus == "rs6967330",]
gsdmb <- nsl[nsl$locus == "rs2305480",]
il33 <- nsl[nsl$locus == "rs928413",]
rad50 <- nsl[nsl$locus == "rs6871536",]
il1r1 <- nsl[nsl$locus == "rs1558641",]



#summary stats below are from the entire datasets (on the server)
summary(cdhr3$percrank)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.8890  0.9529  0.9744  0.9651  0.9876  0.9979
summary(gsdmb$percrank)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.04873 0.16530 0.30740 0.31390 0.43650 0.66590
summary(il33$percrank)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03432 0.61360 0.75440 0.69700 0.91050 0.99690
summary(rad50$percrank)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03091 0.14920 0.26110 0.26160 0.36060 0.60790
summary(il1r1$percrank)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.1115  0.3547  0.4615  0.5031  0.6830  0.8837

asthma.samp <- rbind(cdhr3, gsdmb, il33, rad50, il1r1)
#did this on entire dataset on the server and wrote to file


##################
asthma_nsl <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/asthma_loci_nsl.txt", header=T, stringsAsFactors = TRUE)

asthma_ihs <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/asthma_loci_ihs_26pops.txt", header=T, stringsAsFactors = TRUE)

names(asthma_ihs)[11] <- "ihs_percrank"
names(asthma_nsl)[11] <- "nsl_percrank"

asthma <- inner_join(asthma_ihs, asthma_nsl, by=c("pop", "chr", "locus", "phys_pos", "X1_freq"))


asthma$ihs_percrank <- asthma$ihs_percrank * 100
asthma$nsl_percrank <- asthma$nsl_percrank * 100

asthma$pop <- factor(asthma$pop, levels=rev(unique((asthma$pop)[order(asthma$ord)])))

asthma$snpID <- paste(str_replace(as.character(asthma$chr), "chr", ""), asthma$phys_pos, sep="_")

aaS <- aa[,c("snpID", "REF", "ALT", "derived")]

asthma <- left_join(asthma, aaS, by="snpID")


asthma$DAF <- ifelse(as.character(asthma$derived) == as.character(asthma$ALT), 1-asthma$X1_freq, asthma$X1_freq)
asthma$daf_ihs <- ifelse(as.character(asthma$derived) == as.character(asthma$ALT), -asthma$norm_ihs, asthma$norm_ihs)
asthma$daf_nsl <- ifelse(as.character(asthma$derived) == as.character(asthma$ALT), -asthma$norm_nsl, asthma$norm_nsl)



asthma.nsl.p <- ggplot(asthma, aes(y=pop, x=locus)) + 
  geom_tile(aes(fill=nsl_percrank)) + 
  geom_text(aes(label = paste(round(-daf_nsl, 2), " (", signif(nsl_percrank, 2), "%)", sep="")), size=3) +
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
asthma.nsl.p

asthma.ihs.p <- ggplot(asthma, aes(y=pop, x=locus)) + 
  geom_tile(aes(fill=ihs_percrank)) + 
  geom_text(aes(label = paste(round(-daf_ihs, 2), " (", signif(ihs_percrank, 2), "%)", sep="")), size=3) +
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
asthma.ihs.p

asthma.freq.p <- ggplot(asthma, aes(y=pop, x=locus)) + 
  geom_tile(aes(fill=DAF)) + 
  geom_text(aes(label = round(DAF, 2)), size=3) +
  scale_fill_gradient(low="#EEEEEE", high="#0099CC") +
  #scale_fill_gradient2(low="#0099CC", high = "#33CC33", midpoint=0.5,
  #                     limits=c(0,1)) +
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
asthma.freq.p

#ggsave(plot=asthma.nsl.p, filename = "C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/asthma_loci_nsl_daf.pdf", width=5, height=10, units="in", dpi=300 )

#ggsave(plot=asthma.ihs.p, filename = "C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/asthma_loci_ihs_daf.pdf", width=5, height=10, units="in", dpi=300 )

#ggsave(plot=asthma.freq.p, filename = "C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/asthma_loci_freq_daf.pdf", width=5, height=10, units="in", dpi=300 )


#######################################3
dat <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/samp_selscan_aa_fst.txt", header=F, stringsAsFactors = TRUE)
names(dat) <- c("pop", "chr", "locus", "pos", "pop_freq", "ihs", "ihs_per", "nsl", "nsl_per", "fst_yri", "fst_lwk", "ref", "alt", "aa", "freq")

str(dat)

der.m <- dat[dat$aa == "A" | dat$aa == "T" | dat$aa == "C" | dat$aa == "G",]
der <- filter(dat, aa == "A" | aa == "C" | aa == "G" | aa == "T")

der$pop_daf <- ifelse(as.character(der$aa) == as.character(der$alt), 1-der$pop_freq, der$pop_freq)
der$daf_ihs <- ifelse(as.character(der$aa) == as.character(der$alt), der$ihs, -der$ihs)
der$daf_nsl <- ifelse(as.character(der$aa) == as.character(der$alt), der$nsl, -der$nsl)
der$daf <- ifelse(as.character(der$aa) == as.character(der$alt), 1-der$freq, der$freq)

ggplot(der) + geom_point(aes(x=pop_daf, y=fst_yri, colour=pop))
ggplot(der) + geom_point(aes(x=pop_daf, y=fst_lwk, colour=pop))
ggplot(der) + geom_point(aes(x=pop_daf, y=daf_ihs, colour=pop))
ggplot(der) + geom_point(aes(x=pop_daf, y=daf_nsl, colour=pop))
ggplot(der) + geom_point(aes(x=daf_nsl, y=fst_yri, colour=pop))
ggplot(der) + geom_point(aes(x=daf_nsl, y=fst_lwk, colour=pop))
ggplot(der) + geom_point(aes(x=daf_nsl, y=daf_ihs, colour=pop))

der$snpID <- paste(str_replace(as.character(der$chr), "chr", ""), der$pos, sep="_")
der$snpID <- as.factor(der$snpID)
der$locus <- as.factor(as.character(der$locus))

der2 <- distinct(der)

spread <- der %>%
  group_by(snpID) %>%
  summarise(stdev = sd(pop_daf), avg = mean(pop_daf), freq = mean(daf))

AFvsSD <- ggplot(spread) + 
  geom_point(aes(x=freq, y=stdev)) +
  geom_hline(yintercept=0.08032828, linetype="dashed", colour="red") +
  geom_vline(xintercept=0.815296, linetype="dashed", colour="red")


exp <- left_join(der, spread, by = c("locus"="locus", "daf"="freq"))
ggplot(exp) + geom_point(aes(x=pop_daf, y=stdev))


selspread <- der %>%
  group_by(snpID) %>%
  summarise(avgnsl = mean(daf_nsl), sdnsl = sd(daf_nsl), avgihs = mean(daf_ihs), avgihs = sd(daf_ihs), perN =mean(nsl_per), sdN = sd(nsl_per), perI = mean(ihs_per), sdI = sd(ihs_per))


snpcounts <- der %>% count(snpID)
n26 <- snpcounts[snpcounts$n == 26, 'snpID']

sub <- der %>% filter(snpID %in% n26$snpID)

subspread <- sub %>%
  group_by(snpID) %>%
  summarise(pop_daf=mean(pop_daf, na.rm=TRUE),
            pop_daf_sd=sd(pop_daf, na.rm=TRUE),
            daf_ihs=mean(daf_ihs, na.rm=TRUE),
            daf_ihs_sd=sd(daf_ihs, na.rm=TRUE),
            ihs_per=mean(ihs_per, na.rm=TRUE),
            ihs_per_sd=sd(ihs_per, na.rm=TRUE),
            daf_nsl=mean(daf_nsl, na.rm=TRUE),
            daf_nsl_sd=sd(daf_nsl, na.rm=TRUE),
            nsl_per=mean(nsl_per, na.rm=TRUE),
            nsl_per_sd=sd(nsl_per, na.rm=TRUE),
            fst_yri=mean(fst_yri, na.rm=TRUE),
            fst_yri=sd(fst_yri, na.rm=TRUE),
            fst_lwk=mean(fst_lwk, na.rm=TRUE),
            fst_lwk=sd(fst_lwk, na.rm=TRUE)
            )

subspread <- sub %>%
  group_by(snpID) %>%
  select(snpID, pop_daf, daf_ihs, ihs_per, daf_nsl, nsl_per, fst_yri, fst_lwk) %>%
  summarize_each(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE)))

subspread.long <- melt(subspread, id.vars="snpID")

ggplot(subspread.long, aes(x="stat", y=value)) + 
  facet_wrap(~variable, scales='free_y') + 
  geom_boxplot() +
  geom_point(subspread.long[subspread.long$snpID=="7_105658451",], aes(x="stat", y=value))


ggplot(subspread) + 
  geom_boxplot(aes(x="SD DAF", y=pop_daf_sd)) + 
  geom_hline(yintercept=subspread[subspread$snpID=="7_105658451",'pop_daf_sd']$pop_daf_sd, linetype="dashed", colour="red")
                  
make_boxplot <- function(d, colnum, statstring) {
  p <- ggplot(d) + 
    geom_boxplot(aes(x=statstring, y=stat)) +
    geom_hline(yintercept=d[d$snpID=="7_105658451", statstring]$stat, linetype="dashed", colour="red")
  return(p)
}




ggplot(subspread.long, aes(time,value)) + geom_line() + facet_grid(series ~ .)



daf_nsl_mean.p <- subspread %>%
  gather(-daf_nsl_mean, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = daf_nsl_mean)) +
  geom_point() +
  #stat_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

            

snp.6882 <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/snps_06882_count.txt", header=T, stringsAsFactors = TRUE)

snp.50 <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/snps_050_count.txt", header=T, stringsAsFactors = TRUE)

snp.6882$per <- snp.6882$nn/sum(snp.6882$nn)

n.p <- ggplot(snp.6882, aes(n, per)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Number of Populations") + 
  ylab("Numer of SNPs with DAF >= 0.6882") +
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  #scale_y_continuous(breaks=c(0, 20000, 40000, 60000, 80000, 100000), labels=c("0", "20K", "40K", "60K", "80K", "100K"))
  scale_y_continuous(labels=scales::percent) #+
  #geom_hline(yintercept=0.03732931, linetype="dashed", colour="red")


sumstats <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/snps_06882_sumstats.txt", header=T, stringsAsFactors = TRUE)
names(sumstats)[4] <- "ihs_per_mean"
names(sumstats)[11] <- "ihs_per_sd"


n.ihsper.p <- ggplot(sumstats, aes(x=n, y=ihs_per_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=ihs_per_mean - ihs_per_sd, ymax=ihs_per_mean + ihs_per_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  ylab("iHS percentile") +
  scale_y_continuous(limits = c(0,1), labels=scales::percent) +
  geom_point(aes(x=26, y=mean(der[der$locus=="rs6967330",'ihs_per']), colour="red")) +
  geom_errorbar(aes(x=26, ymin=mean(der[der$locus=="rs6967330",'ihs_per']) - sd(der[der$locus=="rs6967330",'ihs_per']), ymax=mean(der[der$locus=="rs6967330",'ihs_per']) + sd(der[der$locus=="rs6967330",'ihs_per'])), width=.2, colour="red") +
  theme(legend.position="none")

n.nslper.p <- ggplot(sumstats, aes(x=n, y=nsl_per_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=nsl_per_mean - nsl_per_sd, ymax=nsl_per_mean + nsl_per_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  ylab("nSL percentile") +
  scale_y_continuous(limits = c(0,1), labels=scales::percent) +
  geom_point(aes(x=26, y=mean(der[der$locus=="rs6967330",'nsl_per']), colour="red")) +
  geom_errorbar(aes(x=26, ymin=mean(der[der$locus=="rs6967330",'nsl_per']) - sd(der[der$locus=="rs6967330",'nsl_per']), ymax=mean(der[der$locus=="rs6967330",'nsl_per']) + sd(der[der$locus=="rs6967330",'nsl_per'])), width=.2, colour="red") +
  theme(legend.position="none")

n.ihs.p <- ggplot(sumstats, aes(x=n, y=-daf_ihs_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=-daf_ihs_mean - daf_ihs_sd, ymax=daf_ihs_mean + daf_ihs_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  ylab("iHS") +
  geom_point(aes(x=26, y=mean(der[der$locus=="rs6967330",'daf_ihs']), colour="red")) +
  geom_errorbar(aes(x=26, ymin=mean(der[der$locus=="rs6967330",'daf_ihs']) - sd(der[der$locus=="rs6967330",'daf_ihs']), ymax=mean(der[der$locus=="rs6967330",'daf_ihs']) + sd(der[der$locus=="rs6967330",'daf_ihs'])), width=.2, colour="red") +
  theme(legend.position="none")

n.nsl.p <- ggplot(sumstats, aes(x=n, y=-daf_nsl_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=-daf_nsl_mean - daf_nsl_sd, ymax=daf_nsl_mean + daf_nsl_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  ylab("nSL") +
  geom_point(aes(x=26, y=mean(der[der$locus=="rs6967330",'daf_nsl']), colour="red")) +
  geom_errorbar(aes(x=26, ymin=mean(der[der$locus=="rs6967330",'daf_nsl']) - sd(der[der$locus=="rs6967330",'daf_nsl']), ymax=mean(der[der$locus=="rs6967330",'daf_nsl']) + sd(der[der$locus=="rs6967330",'daf_nsl'])), width=.2, colour="red") +
  theme(legend.position="none")

n.popdaf.p <- ggplot(sumstats, aes(x=n, y=pop_daf_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=pop_daf_mean - pop_daf_sd, ymax=pop_daf_mean + pop_daf_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  ylab("DAF") +
  scale_y_continuous(limits = c(0.6,1))

ggplot(sumstats, aes(x=n, y=fst_yri_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=fst_yri_mean - fst_yri_sd, ymax=fst_yri_mean + fst_yri_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), limits=c(0,27)) +
  ylab("Fst YRI") 
  
ggplot(sumstats, aes(x=n, y=fst_lwk_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=fst_lwk_mean - fst_lwk_sd, ymax=fst_lwk_mean + fst_lwk_sd), width=.2) +
  theme_bw() +
  xlab("Number of Populations") + 
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25)) +
  ylab("Fst LWK") 



p1 <- plot_grid(
  n.p,
  n.ihsper.p,
  n.nslper.p,
  ncol = 1,
  align="v"
)

p2 <- plot_grid(
  n.p,
  n.ihs.p,
  n.nsl.p,
  ncol = 1,
  align="v"
)


###########
der <- mutate(der, floor = floor(pop_daf*100))
der <- mutate(der, floor_world = floor(daf*100))


af_pop <- der[der$locus=="rs6967330", c('pop', 'pop_daf', 'floor')]

datalist = list()
for(row in 1:nrow(af_pop)) {
  print(row)
    df <- filter(der, pop==af_pop[row, 'pop'], floor==af_pop[row, 'floor'])
    datalist[[row]] <- df
}

big_data = do.call(rbind, datalist)


df <- read.table("C:/Users/Mary/PepLab/CDHR3/data/asthma_loci/170522_rs6967330AF_stats.txt", header=T, stringsAsFactors = TRUE)

df$col <- ifelse(as.character(df$locus) == "rs6967330", 0, 1)
df$col <- as.factor(df$col)  

p <- ggplot(df, aes(x=pop, y=fst_yri)) +
  geom_boxplot() +
  geom_point(df[df$col==0,], aes(x=pop, y=fst_yri[]))

geom_point( aes(x=factor(cyl)[idx],y=mpg[idx]) )    

ggplot(subset(df, !highlight), aes(x=variable, y=value)) + 
  geom_boxplot() + scale_y_log10() +
  geom_point(                               # add the highlight points
    data=subset(df.2, highlight), 
    aes(x=variable, y=value), 
    color="red", size=5
  )

df$highlight <- ifelse(as.character(df$locus) == "rs6967330", TRUE, FALSE)  

df$pop <- factor(df$pop, levels=rev(unique((asthma$pop)[order(asthma$ord)])))

yri.p <- ggplot(subset(df, !highlight), aes(x=pop, y=fst_yri)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=fst_yri), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("FST - YRI") +
  theme_bw()

lwk.p <- ggplot(subset(df, !highlight), aes(x=pop, y=fst_lwk)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=fst_yri), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("FST - LWK") +
  theme_bw()


daf.p <- ggplot(subset(df, !highlight), aes(x=pop, y=daf)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=daf), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("daf") +
  theme_bw()

nsl.p <- ggplot(subset(df, !highlight), aes(x=pop, y=daf_nsl)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=daf_nsl), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("nSL") +
  theme_bw()

ihs.p <- ggplot(subset(df, !highlight), aes(x=pop, y=daf_ihs)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=daf_ihs), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("iHS") +
  theme_bw()

ihsper.p <- ggplot(subset(df, !highlight), aes(x=pop, y=ihs_er)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=ihs_er), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("% iHS") +
  theme_bw()

nslper.p <- ggplot(subset(df, !highlight), aes(x=pop, y=nsl_per)) + 
  geom_boxplot() + 
  geom_point(                               # add the highlight points
    data=subset(df, highlight), 
    aes(x=pop, y=nsl_per), 
    color="red", size=2
  ) +
  xlab("") +
  ylab("% nSL") +
  theme_bw()

geom_text(data=subset(df, highlight), aes(label=floor))

fst.p <- plot_grid(
  yri.p,
  lwk.p,
  align="v",
  ncol=1
)

sel.p <- plot_grid(
  ihsper.p,
  nslper.p,
  align="v",
  ncol=1
)

fst_sel.p <- plot_grid(
  yri.p,
  lwk.p,
  ihsper.p,
  nslper.p,
  align="v",
  ncol=1
)