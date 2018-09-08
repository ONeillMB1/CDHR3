require(ggplot2)
require(ggtree)
require(phytools)
require(phangorn)
require(ape)
require(Biostrings)
require(seqinr)

cdhr3 <- read.tree("C:/Users/Mary/PepLab/CDHR3/data/RAxML_bestTree.AUTO")
cdhr3.mid <- midpoint(cdhr3)
species.common <- read.tree("C:/Users/Mary/PepLab/CDHR3/data/species_commonname.tree")
species.common.mid <- midpoint(species.common)
species.sci <- read.tree("C:/Users/Mary/PepLab/CDHR3/data/species_scientificname.tree")
species.sci.mid <- midpoint(species.sci)

key <- read.table("C:/Users/Mary/PepLab/CDHR3/data/multiz_100way_key_v2.txt", header=F, sep='\t')


t <- ggtree(cdhr3.mid) + geom_text(aes(label=tip.))

aln <- "C:/Users/Mary/PepLab/CDHR3/data/CDHR3_multiz100way_namesMod.fasta"
aln.sci <- "C:/Users/Mary/PepLab/CDHR3/data/CDHR3_multiz100way_commonnames.fasta"
aln.common <- "C:/Users/Mary/PepLab/CDHR3/data/CDHR3_multiz100way_scinames.fasta"
aln.sci.R <- read.alignment(aln.sci, format="fasta")

species.sci.mid.drop <- drop.tip(species.sci.mid, tip=setdiff(species.sci.mid$tip.label, aln.sci.R$nam))


msaplot(t, aln.sci, window = c(519,539))


t <- ggtree(species.sci.mid.drop) + geom_tiplab(align = TRUE, linetype="dotted")
msaplot(ggtree(species.sci.mid.drop), aln.sci, windo = c(524,534))


species_aln <- msaplot(t, aln.sci, offset = 2, window = c(525,533), color=scale_fill_viridis())
ggsave("C:/Users/Mary/PepLab/CDHR3/data/170411_species_protein_tree_and_alignment.pdf", plot=species_aln, width= 8, height=12)

ggsave(paste("figs/beast/mig_and_bsp/", format(Sys.time(), "%y-%m-%d_"), "full_mig_and_logbsp.pdf", sep = ""), plot=fig2.full, width=9, height=4.5)

species.sci.mid.drop$tip.label <- key[match(species.sci.mid.drop$tip.label, key$V3), 'V7']

cdhr3.mid$tip.label <- key[match(sapply(strsplit(cdhr3.mid$tip.label, "_"), "[[", 3), key$V1), 'V3']

mod.tip.labs <- function(tree, d) {
  tre_tips <- subset(d, isTip == TRUE, c(label, Geography))
  tre_tips$samp <- sapply(strsplit(tre_tips$label, "_"), "[[", 1)
  tre_tips$iso2 <- sapply(strsplit(tre_tips$label, "_"), "[[", 2)
  tre_tips$lin <- sapply(strsplit(tre_tips$label, "_"), "[[", 4)
  tre_tips$newName <- paste(tre_tips$samp, tre_tips$lin, tre_tips$iso2, tre_tips$Geography, sep="_")
  tree@phylo$tip.label <- tre_tips$newName
  return(tree)
}