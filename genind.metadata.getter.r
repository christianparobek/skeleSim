## This is a function to parse a genInd object (from adegenet)
## We want to recover the following information:
##      Number of Pops
##      Number of Loci
##      Number of Samples
##      For each Locus, allele frequency

## Started 18 March 2015
## Started by cp

############################
############################
# 
# if (TRUE){
#   ## Load required packages
#   library(adegenet)
#   
#   ## Make a toy genind object, for haploid SNP data
#   dat <- matrix(sample(c("a","t","g","c"), 25, replace=TRUE), nrow=5)
#   rownames(dat) <- paste("sample.", 1:5)
#   colnames(dat) <- 1:5
#   x <- df2genind(dat, ploidy=1)
#   pop(x) <- c("a","a","b","b","a")
# }

#############################
## Given a genind object...##
#############################

genind.metadata.getter <- function(genind){

## Load required packages ##
stopifnot(require(adegenet))

## Get number of populations
num_pops <- length(genind@pop.names)

## Get number of loci
num_loci <- length(genind@loc.names)

## Get number of samples per population
capture.output(samps_per_pop <- summary(genind)$pop.eff)

## Get ploidy
ploidy <- genind@ploidy

## Get allele freq per locus
freqs <- colMeans(genind@tab)
names(freqs) <- sapply(strsplit(colnames(genind@tab), "\\."), function(x){x[1]})

freq_df <- NULL
freq_df$names <- names(freqs)
freq_df$freqs <- freqs
freq_by_locus <- aggregate(freqs ~ names, freq_df, c)$freqs

## Make a single list to return
genind_metadata <- c(num_pops, num_loci, ploidy)
genind_metadata <- append(genind_metadata, list(as.vector(samps_per_pop)))
genind_metadata <- append(genind_metadata, list(freq_by_locus))

# Name all the elements in our genind_metadata list
names(genind_metadata) <- c("NumberOfPops", "NumberOfLoci", "Ploidy", "SampsPerPop", "FreqByLocus")
# Rename the loci
names(genind_metadata$FreqByLocus) <- 1:num_loci

return(as.list(genind_metadata))
}

