require(dplyr)
require(ggplot2)

dat <- read.table("C:/Users/Mary/Downloads/170421_nSL_superpopOverlap.txt", header=T, sep="\t")

nSL <- dat[,1:6]
rownames(nSL) <- nSL$pos
nSL <- nSL[,2:6]

dat <- tbl_df(dat)

dat <- dat %>%
  rowwise() %>%
  mutate(avg = mean(c(nSLAFR,nSLAMR,nSLEAS,nSLEUR,nSLSAS)), dev = sd(c(nSLAFR,nSLAMR,nSLEAS,nSLEUR,nSLSAS)))


dat <- arrange(dat, dev)

