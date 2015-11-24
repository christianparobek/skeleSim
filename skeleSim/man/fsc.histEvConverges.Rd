\name{fsc.histEvConverges}
\alias{fsc.histEvConverges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fsc.histEvConverges(hist.ev, pop.size, growth.rate, num.mig.mats = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.ev}{
%%     ~~Describe \code{hist.ev} here~~
}
  \item{pop.size}{
%%     ~~Describe \code{pop.size} here~~
}
  \item{growth.rate}{
%%     ~~Describe \code{growth.rate} here~~
}
  \item{num.mig.mats}{
%%     ~~Describe \code{num.mig.mats} here~~
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
function (hist.ev, pop.size, growth.rate, num.mig.mats = NULL) 
{
    if (is.null(hist.ev)) 
        return(TRUE)
    if (length(growth.rate) == 1) 
        growth.rate <- rep(growth.rate, length(pop.size))
    if (length(pop.size) != length(growth.rate)) {
        cat("'pop.size' and 'growth.rate' vectors must be same size.\n")
        return(FALSE)
    }
    if (!fsc.histEvCheck(hist.ev, length(pop.size), num.mig.mats)) 
        return(FALSE)
    hist.ev <- hist.ev[order(hist.ev[, 1], hist.ev[, 2], hist.ev[, 
        3]), , drop = FALSE]
    print(hist.ev)
    for (i in 1:nrow(hist.ev)) {
        gen <- hist.ev[i, 1]
        from <- hist.ev[i, 2] + 1
        to <- hist.ev[i, 3] + 1
        mig <- hist.ev[i, 4]
        new.size <- hist.ev[i, 5]
        pop.size[from] <- pop.size[from] * (1 - mig)
        t <- ifelse(i == 1, gen, gen - hist.ev[i - 1, 1])
        pop.size[to] <- pop.size[to] * exp(growth.rate[to] * 
            t)
        growth.rate[to] <- hist.ev[i, 6]
        pop.size[to] <- pop.size[to] * new.size
        print(pop.size)
    }
    sum(pop.size > 0) == 1
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
