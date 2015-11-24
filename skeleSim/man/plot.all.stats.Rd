\name{plot.all.stats}
\alias{plot.all.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.all.stats(params)
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
    stopifnot(require("reshape2"))
    stopifnot(require("ggplot2"))
    stopifnot(require("gridExtra"))
    results.datafr <- as.data.frame(params@analysis.results)
    num.sc <- length(params@scenarios)
    names.stats <- colnames(params@analysis.results)
    colnames(results.datafr) <- c(names.stats)
    num.stats <- length(names.stats) - 1
    results.melted <- melt(results.datafr, id = "scenario", measure.vars = c(names.stats[-1]))
    ggplot(results.melted, aes(value)) + geom_density() + facet_wrap(scenario ~ 
        variable, ncol = num.stats, scales = "free")
    ggplot(results.melted, aes(value)) + geom_histogram() + facet_grid(variable ~ 
        scenario, scales = "free")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
