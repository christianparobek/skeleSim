\name{rms.run}
\alias{rms.run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rms.run(params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
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
function (params) 
{
    sc <- currentScenario(params)
    skeleland <- rms.init.landscape(num.pops = sc@num.pops, carrying = sc@pop.size, 
        mig.rates = sc@migration, num.loc = sc@num.loci, loc.type = sc@locus.type, 
        mut.rate = sc@mut.rate, seq.length = sc@sequence.length, 
        num.stgs = sc@simulator.params@num.stgs, selfing = sc@simulator.params@selfing, 
        surv.matr = sc@simulator.params@surv.matr, repr.matr = sc@simulator.params@repr.matr, 
        male.matr = sc@simulator.params@male.matr, init.pop.sizes = sc@simulator.params@init.pop.sizes, 
        num.alleles = sc@simulator.params@num.alleles, allele.freqs = sc@simulator.params@allele.freqs)
    if (!is.landscape(skeleland)) 
        stop("landscape not cool")
    skeleland <- landscape.simulate(skeleland, sc@simulator.params@num.gen)
    skeleland_samp <- landscape.sample(skeleland, sc@num.pops, 
        sc@sample.size)
    params@rep.sample <- rms.convert(skeleland_samp, sc@locus.type)
    return(params)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
