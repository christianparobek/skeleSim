\name{landscape.ind.freq}
\alias{landscape.ind.freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
landscape.ind.freq(l, include.states = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{include.states}{
%%     ~~Describe \code{include.states} here~~
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
function (l, include.states = T) 
{
    aml <- vector("list", length(landscape.ploidy(l)))
    for (loc in 1:length(landscape.ploidy(l))) {
        genos <- landscape.locus(loc, l)[, -1:-landscape.democol()]
        ploidy <- landscape.ploidy(l)[loc]
        if ((include.states) & (l$loci[[loc]]$type != 253)) {
            lst <- landscape.locus.states(loc, l)
            names(lst$state) <- lst$aindex
            if (ploidy == 2) {
                genos[, 1] <- unname(lst$state[as.character(genos[, 
                  1])])
                genos[, 2] <- unname(lst$state[as.character(genos[, 
                  2])])
            }
            else {
                genos <- unname(lst$state[as.character(genos)])
            }
        }
        amat <- sapply(as.numeric(names(table(genos))), function(x, 
            genos, pl) {
            if (pl == 2) {
                (genos[, 1] == x) + (genos[, 2] == x)
            }
            else {
                genos == x
            }
        }, genos = genos, pl = ploidy)
        aml[[loc]] <- apply(amat, 2, function(x, pl) {
            x/pl
        }, pl = ploidy)
    }
    do.call(cbind, aml)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
