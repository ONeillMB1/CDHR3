dat <- read.table("C:/Users/Mary/Downloads/HGDPid_populations.txt", header=T, sep='\t')

library(dplyr)

dat <- tbl_df(dat)

pop <- arrange(count(dat, population), desc(n))

reg <- arrange(count(dat, Region), desc(n))

grp <- arrange(count(dat, Pop7Groups), desc(n))
