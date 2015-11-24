\name{rms.convert}
\alias{rms.convert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rms.convert(Rland, locus.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rland}{
%%     ~~Describe \code{Rland} here~~
}
  \item{locus.type}{
%%     ~~Describe \code{locus.type} here~~
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
function (Rland, locus.type) 
{
    ltype <- locus.type
    if (ltype == "microsat") {
        this.rep.result <- landscape.make.genind(Rland)
    }
    else if (ltype == "sequence") {
        states <- as.data.frame(landscape.locus.states(1, Rland))
        genos <- data.frame(pop = landscape.populations(Rland), 
            aindex = Rland$individuals[, 7])
        seq <- merge(genos, states, all.x = T)
        seq <- seq[order(seq$pop), ]
        dna.seq <- strsplit(as.character(tolower(seq$state)), 
            "")
        this.rep.result <- list(strata = data.frame(seq$pop), 
            dna.seq = as.DNAbin(do.call(rbind, strsplit(tolower(as.character(seq$state)), 
                ""))))
    }
    else if (ltype == "SNP") {
        this.rep.result <- landscape.make.genind(landscape.snp.convert(Rland))
    }
    this.rep.result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
