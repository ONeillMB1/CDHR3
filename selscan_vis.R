####iHS
ihs.pop <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/1k_pop_ihs_norm.txt", header=F, sep = "\t")
names(ihs.pop) <- c("pop", "locus_ID", "position", "alt_freq", "ihh1", "ihh0", "iHS", "norm_iHS", "var")

ihs.sup <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/superpop_ihs_norm.txt", header=F, sep = "\t")
names(ihs.sup) <- c("pop", "locus_ID", "position", "alt_freq", "ihh1", "ihh0", "iHS", "norm_iHS", "var")

ihs <- rbind(ihs.sup, ihs.pop)
####nSL
nsl.pop <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/1k_pop_nsl_norm.txt", header=F, sep = "\t")
names(nsl.pop) <- c("pop", "locus_ID", "position", "alt_freq", "sL1", "sL0", "nSL", "norm_nSL", "var")

nsl.sup <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/superpop_nsl_norm.txt", header=F, sep = "\t")
names(nsl.sup) <- c("pop", "locus_ID", "position", "alt_freq", "sL1", "sL0", "nSL", "norm_nSL", "var")

#nsl <- rbind(nsl.sup, nsl.pop) #not enough memory



###Percent rank
require(dplyr)
ihs.tbl <- tbl_df(ihs)

ihs.tbl <- ihs.tbl %>%
  group_by(pop) %>% 
  mutate(percrank=rank(abs(norm_iHS))/length(norm_iHS))

ihs.rs <- filter(ihs.tbl, locus_ID == "rs6967330")

nsl.pop <- tbl_df(nsl.pop)
nsl.sup <- tbl_df(nsl.sup)

nsl.pop <- nsl.pop %>%
  group_by(pop) %>%
  mutate(percrank=rank(abs(norm_nSL))/length(norm_nSL))

nsl.sup <- nsl.sup %>%
  group_by(pop) %>%
  mutate(percrank=rank(abs(norm_nSL))/length(norm_nSL))

nsl.pop.rs <- filter(nsl.pop, locus_ID == "rs6967330")
nsl.sup.rs <- filter(nsl.sup, locus_ID == "rs6967330")

nsl.rs <- rbind(nsl.sup.rs, nsl.pop.rs)

nsl.rs[nsl.rs$pop == "CHB_maf04.nsl.out.100bins.norm", 'pop'] <- "CHB"

#combine datasets

sub.nsl <- select(nsl.rs, pop, norm_nSL, percrank)
sub.ihs <- select(ihs.rs, pop, norm_iHS, percrank)

names(sub.nsl) <- c("pop", "norm_nSL", "nSL_per")
names(sub.ihs) <- c("pop", "norm_iHS", "iHS_per")

sub.nsl$pop <- as.character(sub.nsl$pop)
sub.ihs$pop <- as.character(sub.ihs$pop)

comb <- full_join(sub.ihs, sub.nsl, by="pop")


pops <- read.table("C:/Users/Mary/PepLab/CDHR3/data/selscan_output_summaries/1k_pop.txt", header=F, sep='\t')
names(pops) <- c("pop", "pop_description", "superpop", "ord")

head(comb)
comb2 <- comb %>%
  mutate(pop_des = strsplit(pop, "_")[0])
