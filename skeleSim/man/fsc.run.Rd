\name{fsc.run}
\alias{fsc.run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fsc.run(params)
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
    label <- currentLabel(params)
    sc <- currentScenario(params)
    if (file.exists(params@wd)) 
        for (f in dir(label, full.names = T)) file.remove(f)
    size.mult <- if (sc@locus.type == "sequence") 
        1
    else 2
    file <- fsc.write(num.pops = sc@num.pops, Ne = sc@pop.size * 
        size.mult, sample.size = sc@sample.size * size.mult, 
        sample.times = sc@simulator.params@sample.times, growth.rate = sc@simulator.params@growth.rate, 
        mig.rates = sc@migration, num.chrom = sc@num.loci, hist.ev = as.matrix(sc@simulator.params@hist.ev), 
        locus.params = sc@simulator.params@locus.params, label = label)
    if (.Platform$OS.type == "unix") {
        err <- system(paste(sc@simulator.params@fastsimcoal.exec, 
            "-i", file, "-n 1", ifelse(sc@simulator.params@inf.site.model, 
                "-I", ""), ifelse(params@quiet, "-q", "")), intern = F)
    }
    else {
        err <- shell(paste(sc@simulator.params@fastsimcoal.exec, 
            "-i", file, "-n 1", ifelse(sc@simulator.params@inf.site.model, 
                "-I", ""), ifelse(params@quiet, "-q", "")), intern = F)
    }
    if (err == 0) {
        if (!params@quiet) 
            cat("fastsimcoal exited normally\n")
    }
    else {
        stop("fastsimcoal exited with error ", err, "\n")
    }
    arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
    params@rep.sample <- fsc.read(arp.file, sc)
    params
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
