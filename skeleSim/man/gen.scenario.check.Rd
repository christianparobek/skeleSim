\name{gen.scenario.check}
\alias{gen.scenario.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gen.scenario.check(params)
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
    results.check <- sapply(params@scenarios, function(sc) {
        c(nsizes.eq.npops = length(sc@pop.size) == sc@num.pops, 
            nsamps.eq.npops = length(sc@sample.size) == sc@num.pops, 
            is.mig.square = sapply(sc@migration, function(mig) {
                nrow(mig) == ncol(mig) & nrow(mig) == sc@num.pops
            }), at.lst.1.pop = sc@num.pops >= 1, at.lst.1.loc = sc@num.loci >= 
                1, mut.rate.ok = all((sc@mut.rate >= 0) & (sc@mut.rate < 
                1)), at.lst.1.samp = min(sc@sample.size) > 0, 
            mig.bet.0.1 = sapply(sc@migration, function(mig) {
                all((mig >= 0) & (mig <= 1))
            }), mig.diag.eq.0 = sapply(sc@migration, function(mig) {
                all(diag(mig) == 0)
            }))
    })
    return(results.check)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
