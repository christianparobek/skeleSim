\name{fsc.histEvMat}
\alias{fsc.histEvMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fsc.histEvMat(num.events = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num.events}{
%%     ~~Describe \code{num.events} here~~
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
function (num.events = 0) 
{
    if (num.events == 0) 
        return(NULL)
    hist.ev <- c(num.gen = 0, source.deme = 0, sink.deme = 0, 
        prop.migrants = 1, new.sink.size = 1, new.sink.growth = 0, 
        new.mig.mat = 0)
    do.call(rbind, lapply(1:num.events, function(x) hist.ev))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
