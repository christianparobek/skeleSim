\name{fsc.locusParamsMat}
\alias{fsc.locusParamsMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fsc.locusParamsMat(sc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sc}{
%%     ~~Describe \code{sc} here~~
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
function (sc) 
{
    locus.type <- switch(sc@locus.type, microsat = "MICROSAT", 
        snp = "SNP", sequence = "DNA")
    num.loci <- sc@num.loci
    locus.length <- sc@sequence.length
    mut.rate <- sc@mut.rate
    locus.param.5 <- 1/3
    if (locus.type == "MICROSAT") 
        locus.param.5 <- 0
    if (locus.type == "SNP") 
        locus.param.5 <- NULL
    locus.param.6 <- 0
    if (locus.type \%in\% c("SNP", "DNA")) 
        locus.param.6 <- NULL
    locus.type <- rep(locus.type, num.loci)
    cbind(locus.type = locus.type, locus.length = locus.length, 
        recomb.rate = 0, mut.rate = mut.rate, locus.param.5 = locus.param.5, 
        locus.param.6 = locus.param.6)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
