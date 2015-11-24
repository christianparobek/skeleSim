\name{calc.mratio}
\alias{calc.mratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calc.mratio(gen.data.gtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen.data.gtype}{
%%     ~~Describe \code{gen.data.gtype} here~~
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
function (gen.data.gtype) 
{
    pop.locus.df <- as.matrix(expand.grid(pop = 1:length(unique(gen.data.gtype@strata)), 
        locus = 1:length(gen.data.gtype@loci)))
    calc.mratio.2 <- function(all.freq) {
        if (any(all.freq > 0)) {
            min.all <- min(which(all.freq > 0))
            max.all <- max(which(all.freq > 0))
            sum.all.present <- sum(all.freq > 0)
            sum.all.present/(max.all - min.all + 1)
        }
        else return(NA)
    }
    all_freq_list <- alleleFreqs(gen.data.gtype, by.strata = T)
    mrat_pop_and_loc <- apply(pop.locus.df, 1, function(x) {
        calc.mratio.2(all_freq_list[[x[2]]][, , x[1]][, 1])
    })
    all_freq_list <- alleleFreqs(gen.data.gtype, by.strata = F)
    mrat_loc <- sapply(1:length(all_freq_list), function(x) {
        calc.mratio.2(as.vector(all_freq_list[[x]][, 1]))
    })
    return(c(mrat_loc, mrat_pop_and_loc))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
