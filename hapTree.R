require(ggplot2)
require(ape)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
require("ggtree")
require(phytools)
require(phangorn)

tree <- read.tree("C:/Users/Mary/PepLab/data/CDHR3/smallHap_ind.tree")

mpt <- midpoint(tree)
p <- ggtree(mpt) + geom_tippoint()

p
toCollapse <- which(mpt$node.label == "")
nodes <- as.numeric(mpt$node.label)
toCollapse2 <- which(is.na(nodes))
toCollapse3 <- which(nodes < 0.5)


for (i in toCollapse2) {
  p <- p %>% collapse(node = i)
  return(p)
}


#ggtree(mpt) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

ggtree(mpt)

cp <- ggtree(mpt) %>% collapse(node=5825)
new <- cp + geom_point2(aes(subset=(node == 5825)), size=5, shape=23, fill="steelblue")
new %>% collapse(node=1056)
new <- new + geom_point2(aes(subset=(node == 1056)), size=5, shape=23, fill="steelblue")
