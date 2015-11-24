\name{fsc.histEvCheck}
\alias{fsc.histEvCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fsc.histEvCheck(hist.ev, num.pops, num.mig.mats = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist.ev}{
%%     ~~Describe \code{hist.ev} here~~
}
  \item{num.pops}{
%%     ~~Describe \code{num.pops} here~~
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
function (hist.ev, num.pops, num.mig.mats = NULL) 
{
    if ((!is.numeric(as.matrix(hist.ev))) & (!(is.matrix(hist.ev) | 
        is.data.frame(hist.ev)))) {
        cat("'hist.ev' must be a numerical matrix or numerical dataframe.\n")
        return(FALSE)
    }
    if (ncol(hist.ev) != 7) {
        cat("'hist.ev' must have 7 columns.\n")
        return(FALSE)
    }
    if (any(hist.ev[, 2] >= num.pops | hist.ev[, 3] >= num.pops)) {
        cat("values in columns 2 and 3 (source and sink demes) cannot be greater than the number of populations.\n")
        return(FALSE)
    }
    if (any(hist.ev[, 4] < 0 | hist.ev[, 5] < 0)) {
        cat("values in columns 4 (proportion of migrants) and 5 (new size for sink deme) must be greater than 0.\n")
        return(FALSE)
    }
    if (!is.null(num.mig.mats)) 
        if (any(hist.ev[, 6] > num.mig.mats)) {
            cat("values in column 6 (new migration matrix) cannot be greater than the number of migration matrices.\n")
            return(FALSE)
        }
    TRUE
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
