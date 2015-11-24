\name{overall_stats}
\alias{overall_stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overall_stats(results_gtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results_gtype}{
%%     ~~Describe \code{results_gtype} here~~
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
function (results_gtype) 
{
    ovl <- overallTest(results_gtype, nrep = 5, quietly = TRUE)
    ovl.result <- ovl$result[complete.cases(ovl$result[, 1]), 
        ]
    pnam <- c()
    for (i in 1:nrow(ovl.result)) pnam <- c(pnam, rownames(ovl.result)[i], 
        paste(rownames(ovl.result)[i], "pval", sep = ""))
    global.wide <- as.vector(t(ovl.result))
    names(global.wide) <- pnam
    global.wide
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
