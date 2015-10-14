####ANALYSIS

#########################################################################################
##HAPLOID SEQUENCE DATA
#########################################################################################
require(ape)

# require(adegenet)
# adegenet from Thibaut's github
  # library(devtools)
  #install_github("thibautjombart/adegenet")
library("adegenet")

# Eric's strataGdevel
  # library(devtools)
  # install_github("ericarcher/swfscMisc/swfscMisc")
  # install_github("ericarcher/strataG.devel/strataG.devel")
library(strataGdevel)

# require(strataG)
library(poppr)
library(MASS)

# Need Thibaut's apex from repository
  # library(devtools)
  # install_github("thibautjombart/apex")
library(apex)

##Empirical and simulated data ('x') should be in DNAbin format, and a metadata df should
##also be created (i.e. 'strata' three columns: id.col(indivID), strata.col(pops), locus.col(haps))
##depending on required summary stats, these files will be coverted to either strataG gtype and/or
###adegenet genind and then to hierfstat format

##### Testing with strataG.devel
# library(strataG.devel)

# Do we need an empty params object?
#params <- list()

#################################  format rep.results from multidna to gtypes #########
# load test multidna.rdata into global environment
class(rep.result)

# multidna, from apex
class(rep.result[[2]])

genes <- rep.result[[2]]#the multidna object
names(genes@dna) <- paste("gene", 1:length(genes@dna))
id <- genes@labels
df <- data.frame(id = id, strata = rep.result[[1]])
gene.labels <- matrix(id, nrow = length(id), ncol = getNumLoci(rep.result[[2]]))
colnames(gene.labels) <- paste("gene", 1:ncol(gene.labels), sep = "_")
df <- cbind(df, gene.labels)
test.g <- df2gtypes(df, 1, sequences = genes)
summary(test.g)

nLoc(test.g)
locNames(test.g)

length(genes)
rep.result
length(rep.result$dna.seqs@dna)

# where does this come from?
#length(seq.names)

# Where do I get params from?
#params@analysis.result <- rep.result
num_loci <- getNumLoci(rep.result[[2]])


## will get a gtypes object
data(dolph.msats)
data(dolph.strata)
msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, all.y = TRUE)
msats <- df2gtypes(msats.merge, ploidy = 2)

# msats <- dolph.msats
alleleFreqs <- alleleFreqs(msats, by.strata = TRUE)
by.loc <- sapply(alleleFreqs, function(loc) {
  mat <- loc[, "freq", ]
  rowSums(apply(mat, 1, function(r) {
    result <- rep(FALSE, length(r))
    if(sum(r > 0) == 1) result[r > 0] <- TRUE
    result
  }))
})
rownames(by.loc) <- strataNames(msats)
perLocus <- colSums(by.loc) #this has the number of alleles that are private per locus
t(by.loc) #the rows will be have the private alleles for each population by locus




#setwd("C:/Users/deprengm/Dropbox/Hackathon/Datasets")
#df <- read.fasta("name_checked_acanigror_CyB_JD.fasta")
#str(df)

class(df)

#if(class(df) == "DNAbin"){
#  df.genind <- multidna2genind(df)
#  class(df.genind)
#  df.gt <- genind2gtypes(df)
#}

### July 13, 2015  ###########################################################
### sequence data, vector of pop assignments and DNAbin objects,
#   NO!: a list with population column and sequence ("hap") column data.frame, column name "strata" (one column)
#   value is a label associated with sequence in DNAbin object
### Yes: need to colapse sequence into haplotype, need vector
### Yes: list one element population assigments, other vector of sequences,
##        list(strata = pop.data[,1], dna.seq = as.DNAbin(dna.seq)) <- that's an apex, not DNAbin which is a class of a list of DNAbins

## then don't need to worry about genind objects, no conversions

#
# sequence data
# step one, convert to DNAbin? take DNAbin?

# data out of simulation will be the list,
####  1. vector of population assigments

# use pairwiseTest, strip character vectors, just numbers, add in the overall
# then know replicates for z dimension
# set array from number of populations, if pairwise + 1 for overall, and then
#   runs for simulation   <- from load params when you tell it number of pops



pairwiseTest.out <- pairwiseTest(msats)
pairwiseTest.out$result

smsat <- summary(msats)
str(smsat)

# strataG uses this order to pick names, can use this to generate
# population 1v1, 1v2 ...
?combn
combn(1:3, 2)
combn(1:9, 2)

## each run of the simulation would produce a new matrix from
# such analayses as pairwiseTest to array
# No: xyzArray <- sapply(1:nrep, pairwiseTest.out, simplify = 'array')

foo <- combn(1:5, 2)
names <- paste(foo[,1], collapse = "v")

rnames <- c("overall", apply(foo, 2, function(x) paste(x, collapse = "v")))
rnames

#####################################################################
####      ROW NAMES
# for pairwise populations, take number of populations = 'npop'
npp <- combn(1:npop, 2)
names <- c("overall", apply(npp, 2, function(x) paste(x, collapse = "v")))

# for pairwise loci, take number of loci = nloc
npl <- combn(1:nloc, 2)
names <- c("overall", apply(npl, 2, function(x) paste(x, collapse = "v")))

# for simple populations
names <-  c("overall", 1:npop)

# for loci
names <- c("overall", 1:nloc)
#####################################################################
## to install from github: install_github("jgx65/hierfstat")  Jerome's

# testing ones
analyses <- c("allel","Freq","prop")
nrep <- 5


params
#params@num.reps <- 1

## Make the array from the load params info
# nanalyze <- the number of analyses we'll do
npop <- 3 # npop <- the number of rows, pops or pairwise or loci plus overall
# nrep <- the number of replicates
# dimension names will come from the choice of pairwise, by
#   by population, by loci
population.by.loci <- array(0, dim = c(length(names), length(analyses), params@num.reps),
                            dimnames = list(names, analyses, 1:nrep))

## array for global - summary of loci for each population
overall.population <- array(0, dim = (length(names), length(analyses), params@num.reps))

## array for

## test making list of $Global, $Pops, $locus, $pws which each could hold list
# per scenario of analysis columns and sample rows

an.req <- c(Global = T, Pops = F, locus = T, pws = F)

for(x in names(which(an.req)))


####################################################################
# now take each run with all the analyses and put into one matrix (per replicate)

# first do multiDNA to gtypes and genind to gtypes
# params@analysis.result will be multidna
# example rep.result
rep.result
class(rep.result)
# A multidna with two genes
length(rep.result$dna.seqs@dna)

data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x

length(x@dna)

# example dolphin
msats
str(locNames(msats))


str(msats)
length(msats@loci)

#Subset is removed from strataG, instead use gtype[id, locus, strata] to subset
#lapply(locNames(msats), function(x){
#  gtypes_1 <- subset(msats, loci = x)
#  ovl <- overallTest(gtypes_1, nrep = 5, quietly = TRUE)
#})


########### HERE ######### Need to collapse/concatenate lists
ovl.loc <- lapply(locNames(msats), function(x){
  gtypes_1 <- msats[,x,]  #[individuals, loci, strata]
  ovl <- overallTest(gtypes_1, nrep = 5, quietly = TRUE)
  ovl$result
})


ovl.out <- do.call(cbind, ovl.loc)

pnam <- c()
for(i in 1:length(colnames(global))){
  pnam <- c(pnam,paste(colnames(global)[i],"pval", sep = ""))
}
global.wide <- c(global[1,],global[2,])
names(global.wide) <- c(colnames(global),pnam)



overall_stats(msats)
################################
# for loop with returning container - the last line of the function, lapply always returns a list, sapply will simplify by end value

mat <- t(sapply(locNames(msats), function (l){
  gtypes_this_loc<-msats[,l,]
  overall_stats(gtypes_this_loc)
}))
###############################

scenario.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)

#num_loci is known from the parameters gathered earlier

num_loci <- nLoc(msats)
num_reps <- 5
analyses <- names(mat[2,])
num_analyses <- length(analyses)

scenario.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci,num_analyses,num_reps),
                                               dimnames = list(1:num_loci,analyses,1:num_reps))

scenario.results[[group]][[curr_scn]] <- array(0, dim=c(1+15,14,100))


which(locNames(msats) == locNames(msats)[2])
#mat[2,] <- 1:14

# example to test
genes <- rep.result[[2]]
#genes is a list
class(genes)
names(genes@dna) <- paste("gene", 1:length(genes@dna))
id <- genes@labels
df <- data.frame(id = id, strata = rep.result[[1]], hap = id)
class(df)
# errors now with designating sequences = genes
test.g <- df2gtypes(df, 1) #, sequences = genes
class(test.g) # multiDNA to gtypes


# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1894623/
# need apex multiDNA object

### Global, Genotypes
ovl <- overallTest(msats, nrep = 5, stat.list = statList("chi2"), quietly = TRUE)
t(ovl$result)
# need new columns with p-value for each
global <- t(ovl$result)
as.vector(t(ovl$result)) # by row to a vector

# Eric's example
data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x.g <- sequence2gtypes(x)
strata(x.g) <- c("A", "B")

# for a multiDNA object, need to add row for results of second gene
  ovl.multi <- sapply(locNames(x.g), function(n) {
    ovl.raw <- overallTest(x.g[,n,], nrep=5, stat.list=statList("chi2"), quietly = TRUE)
    ovl.raw$result
  })

t(ovl.multi)

#
ovl.multi <- sapply(locNames(x.g), function(n) {
  ovl.raw <- overallTest(x.g[,n,], nrep=5, stat.list=statList("chi2"), quietly = TRUE)
  ovl.raw$result
})



# vs.
ovl.2 <- overallTest(x.g[,"gene1",], nrep = 5, stat.list = statList("chi2"), quietly = TRUE)
t(ovl.2$result)
# need new columns with p-value for each
global <- t(ovl.2$result)


rownames(ovl.multi)
rownames(ovl.2$result)

ovl.2$result[which(rownames(ovl.2$result),"Chi2"),]

#
pnam <- c()
for(i in 1:length(rownames(ovl$result))){
    pnam <- c(pnam,rownames(ovl$result)[i],paste(rownames(ovl$result)[i],"pval", sep = ""))
}
#global.wide <- c(ovl$result[1,],ovl$result[2,])
global.wide <- c(global[1,],global[2,])
names(global.wide) <- c(rownames(ovl$result),pnam)

ovl <- overallTest(test.g, nrep = 5, stat.list = statList("chi2"), quetly = TRUE)
t(ovl$result)
# need new columns with p-value for each
global <- t(ovl$result)
pnam <- c()
for(i in 1:length(colnames(global))){
  pnam <- c(pnam,paste(colnames(global)[i],"pval", sep = ""))
}
global.wide <- c(global[1,],global[2,])
names(global.wide) <- c(colnames(global),pnam)


# Population, Genotypes
alfre <- alleleFreqs(msats, by.strata = TRUE)
# use rep.result to test
alfre <- alleleFreqs(test.g)  # can skip this
# may want an allele frequency spectrum, distirbution of these
#   made up 'neutral' alleles
#   if we do a 'pre summary' we lose this information -
#   option to write your own, you can save data
#   if rep.result is a list, anyone can add another
#     automate it so user can add any other
# list of frequency of allels per locus per population
# might take a long time, do we need that many columns?
str(alfre)
alfre$hap
length(alfre[[1]][,,1])
alfre[[1]][,,1] 3 first DNAbin in multidna freq and prop of alleles per locus
for(i in 1:length(alfre))

## remove chi2 but keep the p value
#  no.alleles, allelic.richness, obs.He, exp.He, theta, unique.alleles
smry <- summarizeLoci(msats, by.strata = TRUE)


library(reshape2)
melt.smry <- melt(smry[[1]])

smry.wide <- t(melt.smry[,3])
colnames(smry.wide) <- paste(melt.smry[,1],melt.smry[,2],sep="")


smry <- summarizeLoci(test.g, by.strata = TRUE)
colnames(summarizeLoci(msats))

########################################################################################
########################################################################################
########################################################################################
## Need to add if genes > 1 then {
# Pull analyses names from ...
smry.multi <- t(sapply(locNames(msats), function(n) {
  unname(summarizeLoci(msats[, n, ]))
}))

# extract column names and reformat multigene... wouldn't have without multiple loci...
multi.locus <- sapply(locNames(x.g), function(n) {
  summarizeLoci(x.g[,n,])
})

l.nam <- locNames(x.g)
c.nam <- colnames(summarizeLoci(x.g[,l.nam[1],]))

# formatting works but doesn't make sense to summarizeLoci over a sequence....
multi.t.locus <- t(multi.locus)
dimnames(multi.t.locus)[[2]] <- c.nam

# Can't remember now, I want a row per individual and a long vector of resutls... right?
dimnames(smry.multi)[[2]] <- colnames(summarizeLoci(msats))

# Why unname??
n <- locNames(msats)[1]
smry.multi <- sapply(locNames(msats), function(n) {
  summarizeLoci(msats[, n, ])
})



####
melt.smry <- melt(smry[[1]])
smry.wide <- t(melt.smry[,3])
colnames(smry.wide) <- paste(melt.smry[,1],melt.smry[,2],sep="")


### global, haplotype
ovl <- overall



### pairwise, haplotypes
# test data are: test.g and smsat
# x.g from multigene sapply example
# want as wide as there are alleles
# fairly slow:
pws <- pairwiseTest(msats, nrep = 5, stat.list = list(statGst, quietly = TRUE))
pws
test.y <- pws[[1]][-c(1:5)]
pws[[2]][-c(1:5)]
nrow(test.y)
length(pws[[1]][1])
rownames(test.y) <- as.matrix(pws[[1]][1])

pws.out <- pws$result[-c(1:5)]
rownames(pws.out) <- as.matrix(pws[[1]][1])
pws.out

#multigene
pws.multi <- sapply(locNames(x.g), function(n) {
  pairwiseTest(x.g[, n, ], nrep = 5, stat.list = list(statGst, quietly = TRUE))
})

pws.multi.out <- pws.multi[[1]][-c(1:5)]
rownames(pws.multi.out) <- as.matrix(pws.multi[[1]][1])

# What happens to an example with one gene
data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
wood.g <- sequence2gtypes(x)
strata(wood.g) <- sample(c("A", "B"), nInd(wood.g), rep = TRUE)
wood.g
gene1 <- wood.g[, "gene1", ]
gene1.dnabin <- getSequences(sequences(gene1))
class(gene1.dnabin)



pws.multi <- sapply(locNames(msats), function(n) {
  pairwiseTest(x.g[, n, ], nrep = 5, stat.list = list(statGst, quietly = TRUE))
})

pws.multi.out <- pws.multi[[1]][-c(1:5)]
rownames(pws.multi.out) <- as.matrix(pws.multi[[1]][1])


#
sA <- sharedAlleles(msats)[,-c(1:2)]
sA <- sharedAlleles()
names(sA)
nsharedAlleles <- paste("sharedAlleles", names(sA), sep = ".")
names(sA) <- nsharedAlleles

pws.out <- cbind(pws.out, sA)

sharedAlleles(msats)
sharedAlleles(test.g)
# http://www.genetics.org/content/197/2/769.full.pdf
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3352520/
# http://www.genetics.org/content/197/2/769.full.pdf
# http://mbe.oxfordjournals.org/content/early/2012/10/03/molbev.mss232.full.pdf
# http://www.molecularecologist.com/2015/04/visualizing-linkage-disequilibrium-in-r/
# http://stats.stackexchange.com/questions/23117/which-mean-to-use-and-when
# emperical data alongside and want to us abc simulation to make sure everything calculated the same way

class(test.g)

### genotype is genind: chord distance,


### haplotype is class(multiDNA) class the object, Nei'sDa (distance)
### params has a data type skelesim
## genind genotype
### need to work over all DNAbin in multiDNA: multidna@dna which is a list of DNAbin matrices


# gtypes object already has per population:
#       N, average num alleles, percent unique, heterozygosity
#                           per locus:
#       number genotyped, num allels, percent unique, obs heterozygosity

# Some data about populations
summary(msats)$strata.smry

# Shared alleles
# strataG
# ? how to turn this into one value per population??
sharedAlleles(msats) # per locus
propSharedLoci(msats) # per population columns 3:5



# HWE
# HWEgenepop(msats) ## required GENEPop
# library(pegas)
genindobj <- gtypes2genind(msats)
hw.test(genindobj) # B = # replicates for MC or regular test B = 0
hw.test(genindobj, B = 0)
# library(genetics)
# HWE.test(genotype) # some other data type...



# Some data about loci
summary(msats)$locus.smry


#####################################################################

#foo.array <- array(0, c(2,3,4))

# Then fill it in as the analyses are being done
#foo.array[,,1] <- matrix(1:3, nrow = 2, ncol = 3)

for(i in 1:nrep){
  sim.data.analysis[,,i] <- matrix(all the data)
}

######################################################################







###################################################### end July 13, 2015




















##### Start March Durham  ############################################
####    Old!! From Durham, will now use mostly strataG
#FILE CONVERSION
##Convert to genind, then hierfstat

# From adegenet
# xgenind <- DNAbin2genind(df, pop=strata$pops, exp.char=c("a","t","g","c"))
xgenind <- DNAbin2genind(df, exp.char=c("a","t","g","c"))

xhierfstat <- genind2hierfstat(xgenind, pop=strata$pops)

##Convert to gtype - convert DNAbin file (x) to x.haps, strata will be a df in a list with another df of dna.seq
DNAbin2gtypes <- function(go, strata, popcol = 1){

  if(length(stratafoo > 1)){
    sv <- which(names(strata) == popcol)
  } else {
    sv <- 1
  }

  go <- as.matrix(go)

  x <- sapply(rownames(go), function(n) {
    as.character(go[n, ])[1, ]
  }, simplify = FALSE)

  gtypes(x, strata.vec = sv)

}


#SUMMARY STATISTICS
##population specific
## from genind2hierfstat above
localFst <- betai(xhierfstat)$betaio  #local Fst (Beta of Weir and Hill 2002 NOT of Foll and Gaggiotti 2006)


#### HapFreq <- hap.freqs(gtypes.out) #haplotype frequencies of each population
### Once they all work, make one function that allows choice for which to
##  run, also choice if by population or not - won't that be faster
##  if there isn't population sturcture to, to do by whole group
HapDiv <- haplotypic.diversity(gtypes.out) #haplotype diversities of each population
NucDiv <- nucleotide.diversity(gtypes.out, bases = c("a", "c", "g", "t")) #nucleotide diversity by site
UHap <- pct.unique.haplotypes(gtypes.out)
TajD <- tajimas.d(gtypes.out)

##population pair-wise
NeiD <- nei.Da(gtypes.out) #Nei's Da for each pair of populations
Fst <- pairwise.test(gtypes, stats = "fst", nrep = 10000, keep.null = FALSE,
                     num.cores = 1, quietly = FALSE, write.output = FALSE, ...)
Phist <- pairwise.test(gtypes, stats = "phist", nrep = 10000, keep.null = FALSE,
                       num.cores = 1, quietly = FALSE, write.output = FALSE, ...) #Excoffier et al. 1992

#########################################################################################
##MSAT/SNP DATA
#########################################################################################
require(adegenet)
require(mmod)
require(hierfstat)

##Empirical and simulated data ('x') should be in adegenet genind format, and a metadata df should also be created (i.e. 'strata' three columns: id.col(indivID), strata.col(pops), locus.col(haps))
##depending on required summary stats, these files will be coverted to either genpop (which can be used by adegenet and mmod) or hierfstat data format

#convert to genpop
genind2genpop(x,pop=strata$pops,missing=c("NA"),quiet=FALSE,process.other=FALSE, other.action=mean)

#convert to hierfstat format
genind2hierfstat(x,pop=strata$pops)


#SUMMARY STATISTICS
##population specific
localFst <- betai(xhierfstat) #local Fst (Beta of Weir and Hill 2002 NOT of Foll and Gaggiotti 2006)
HWE <- HWE.test(xgenind,pop=strata$pops,permut=10000,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix")) #DECIDE WHETHER full OR matrix is better output, ALSO HOW MANY permutation and SIMULATIONS
AlFreqs <- scaleGen(xgenpop, center=TRUE, scale=TRUE,
                    missing=c("NA","0","mean"),truenames=TRUE) #computes scaled allele frequencies
AlRich <- allelic.richness(xhierfstat,min.n=NULL,diploid=TRUE) #rarified allele counts, produces a table of rows (loci) and columns (pops), min.all is the number of alleles used for rarefaction

##population pair-wise
DJost <- pairwise_D(xgenind, linearized = FALSE) #This function calculates Jost's D
GprimeST <- pairwise_Gst_Hedrick(xgenind, linearized = FALSE) #This function calculates Hedrick's G'st
GST <- pairwise_Gst_Nei(xgenind, linearized = FALSE) #This function calculates Nei's Gst
FST <- pairwise.fst(xgenind, pop=strata$pops, res.type=c("dist","matrix"), truenames=TRUE) #Nei's Fst (Nei 1973) Ht - ((Hs(A) + Hs(B)/(n_A+n_B)) / Ht ) DECIDE WHETHER dist OR matrix is better output
