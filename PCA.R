#install.packages("irlba")
#install.packages("threejs")

require(irlba)
require(threejs)

setwd("C:/Users/Mary/PepLab/data/CDHR3/PCA/")


#p = pipe("zcat CDHR3.vcf.gz | sed /^#/d | cut -f '10-' | ./a.out | cut -f '1-2'")
x = read.table("p", colClasses=c("integer","integer"), fill=TRUE, row.names=NULL)

# Convert to a sparse matrix of people (rows) x variant (columns)
cdhr3 = sparseMatrix(i=x[,2], j=x[,1], x=1.0)

# Inspect the dimensions of this matrix
print(dim(cdhr3))

cm = colMeans(cdhr3)
p = irlba(cdhr3, nv=3, nu=3, tol=0.1)

#scatterplot3js(p$u)

require(scatterplot3d)
scatterplot3d(p$u)

# Read just the header of the chromosome file to obtain the sample identifiers
ids = readLines("ids")

# Download and parse the superpopulation data for each sample, order by ids
ped = read.table("integrated_call_samples_v2.20130502.ALL.ped",sep="\t",header=TRUE,row.names=2)[ids,6,drop=FALSE]

# Download the subpopulation and superpopulation codes
# WARNING: These links occasionally change. Beware!
#pop = read.table("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv",sep="\t",header=TRUE)
#pop = pop[1:26,]
#super = pop[,3]
#names(super) = pop[,2]
#super = factor(super)
# The last rows of pop are summary data or non-relevant:

# Map sample sub-populations to super-populations
#ped$Superpopulation = super[as.character(ped$Population)]

# Plot with colors corresponding to super populations
N = length(levels(ped$Population))
scatterplot3d(p$u, color=rainbow(N)[ped$Population], size=0.5)
