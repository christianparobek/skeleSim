\name{rms.init.landscape}
\alias{rms.init.landscape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rms.init.landscape(num.pops = NULL, carrying = NULL, sample.size = NULL, mig.rates = NULL, num.loc = NULL, loc.type = NULL, mut.rate = NULL, seq.length = NULL, num.stgs = NULL, selfing = NULL, surv.matr = NULL, repr.matr = NULL, male.matr = NULL, init.pop.sizes = NULL, num.gen = NULL, num.alleles = NULL, allele.freqs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num.pops}{
%%     ~~Describe \code{num.pops} here~~
}
  \item{carrying}{
%%     ~~Describe \code{carrying} here~~
}
  \item{sample.size}{
%%     ~~Describe \code{sample.size} here~~
}
  \item{mig.rates}{
%%     ~~Describe \code{mig.rates} here~~
}
  \item{num.loc}{
%%     ~~Describe \code{num.loc} here~~
}
  \item{loc.type}{
%%     ~~Describe \code{loc.type} here~~
}
  \item{mut.rate}{
%%     ~~Describe \code{mut.rate} here~~
}
  \item{seq.length}{
%%     ~~Describe \code{seq.length} here~~
}
  \item{num.stgs}{
%%     ~~Describe \code{num.stgs} here~~
}
  \item{selfing}{
%%     ~~Describe \code{selfing} here~~
}
  \item{surv.matr}{
%%     ~~Describe \code{surv.matr} here~~
}
  \item{repr.matr}{
%%     ~~Describe \code{repr.matr} here~~
}
  \item{male.matr}{
%%     ~~Describe \code{male.matr} here~~
}
  \item{init.pop.sizes}{
%%     ~~Describe \code{init.pop.sizes} here~~
}
  \item{num.gen}{
%%     ~~Describe \code{num.gen} here~~
}
  \item{num.alleles}{
%%     ~~Describe \code{num.alleles} here~~
}
  \item{allele.freqs}{
%%     ~~Describe \code{allele.freqs} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (num.pops = NULL, carrying = NULL, sample.size = NULL, 
    mig.rates = NULL, num.loc = NULL, loc.type = NULL, mut.rate = NULL, 
    seq.length = NULL, num.stgs = NULL, selfing = NULL, surv.matr = NULL, 
    repr.matr = NULL, male.matr = NULL, init.pop.sizes = NULL, 
    num.gen = NULL, num.alleles = NULL, allele.freqs = NULL) 
{
    if (exists("skeletonland")) 
        rm(skeletonland)
    skeletonland <- landscape.new.empty()
    skeletonland <- landscape.new.floatparam(skeletonland, s = selfing)
    skeletonland <- landscape.new.intparam(skeletonland, h = num.pops, 
        s = num.stgs, cg = 0, ce = 0, totgen = 1e+05, maxland = 1e+05)
    skeletonland <- landscape.new.switchparam(skeletonland, re = 0, 
        rd = 0, mp = 1, dd = 0)
    for (i in 1:num.pops) skeletonland <- landscape.new.local.demo(skeletonland, 
        surv.matr, repr.matr, male.matr)
    epoch_s_matr <- matrix(0, nrow = 4, ncol = 4)
    epoch_r_matr <- landscape.mig.matrix(h = num.pops, s = num.stgs, 
        mig.model = "stepping.stone.linear")$R
    epoch_m_matr <- matrix(0, nrow = 4, ncol = 4)
    skeletonland <- landscape.new.epoch(skeletonland, epochprob = 1, 
        epoch_s_matr, epoch_r_matr, epoch_m_matr, startgen = 0, 
        extinct = NULL, carry = carrying)
    rms.locus.type = NULL
    if (loc.type == "SNP") {
        rms.locus.type = 2
        num.alleles = 4
        seq.length = rep(1, num.loc)
    }
    if (loc.type == "microsat") 
        rms.locus.type = 1
    if (loc.type == "sequence") 
        rms.locus.type = 2
    for (l in 1:num.loc) skeletonland <- landscape.new.locus(skeletonland, 
        type = 1, ploidy = 2, mutationrate = mut.rate[l], transmission = 0, 
        numalleles = num.alleles[l], frequencies = allele.freqs, 
        allelesize = seq.length[l])
    skeletonland <- landscape.new.individuals(skeletonland, init.pop.sizes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
