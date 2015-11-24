\name{summ.stats.table}
\alias{summ.stats.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summ.stats.table(params)
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
    results.datafr <- params@analysis.results
    num.sc <- length(params@scenarios)
    names.stats <- colnames(params@analysis.results)
    colnames(results.datafr) <- c(names.stats)
    num.stats <- length(names.stats) - 1
    table.means <- data.frame(matrix(NA, nrow = num.sc, ncol = num.stats))
    table.sd <- data.frame(matrix(NA, nrow = num.sc, ncol = num.stats))
    colnames(table.means) <- names.stats[-1]
    colnames(table.sd) <- names.stats[-1]
    for (i in 1:num.stats) table.means[, i] <- tapply(results.datafr[, 
        1 + i], results.datafr[, 1], mean)
    for (i in 1:num.stats) table.sd[, i] <- tapply(results.datafr[, 
        1 + i], results.datafr[, 1], sd)
    results.list <- list()
    results.list$means <- table.means
    results.list$sd <- table.sd
    params@summary.results <- results.list
    return(params)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
