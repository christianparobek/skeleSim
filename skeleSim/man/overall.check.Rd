\name{overall.check}
\alias{overall.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overall.check(params)
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
    params@other.checks <- non.scenario.check(params)
    print(params@other.checks)
    prv_chk <- params@sim.scen.checks
    ths_chk <- rbind(params@sim.check.func(params), gen.scenario.check(params))
    print(prv_chk)
    print(ths_chk)
    if (is.null(prv_chk)) 
        params@sim.scen.checks <- ths_chk
    else {
        for (i in rownames(ths_chk)) {
            if (i \%in\% rownames(prv_chk)) 
                prv_chk[i, ] <- ths_chk[i, ]
            else {
                rbind(prv_chk, ths_chk[i, ])
                rownames(prv_chk)[nrow(prv_chk)] <- i
            }
        }
        params@sim.scen.checks <- prv_chk
    }
    return(params)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
