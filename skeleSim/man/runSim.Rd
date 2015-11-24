\name{runSim}
\alias{runSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runSim(params, num.secs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{num.secs}{
%%     ~~Describe \code{num.secs} here~~
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
function (params, num.secs = NULL) 
{
    params <- overall.check(params)
    if (!all(params@other.checks, params@sim.scen.checks)) {
        filewr <- "error.log"
        write.csv(params@sim.scen.checks, file = filewr)
        print(params@other.checks)
        write("\n", file = filewr, append = T)
        write.table(params@other.checks, file = filewr, append = T)
        stop("parameters do not pass checks; see error log for details")
    }
    print("params checked")
    if (file.exists(params@wd)) {
        unlink(params@wd, recursive = TRUE, force = TRUE)
    }
    dir.create(params@wd)
    wd <- setwd(params@wd)
    results <- list(timing = list())
    params@analysis.results <- NULL
    tryCatch({
        num.sc <- length(params@scenarios)
        num.reps <- params@num.reps
        sc.vec <- rep(1:num.sc, num.reps)
        rep.vec <- rep(1:num.reps, each = num.sc)
        params@scenario.reps <- cbind(scenario = sc.vec, replicate = rep.vec)
        quit <- FALSE
        num.iter <- nrow(params@scenario.reps)
        results$timing$start.time <- Sys.time()
        for (i in 1:num.iter) {
            if (quit) 
                break
            params@current.scenario <- params@scenario.reps[i, 
                "scenario"]
            params@current.replicate <- params@scenario.reps[i, 
                "replicate"]
            params <- params@sim.func(params)
            params <- params@rep.analysis.func(params)
            label <- currentLabel(params)
            file <- paste(label, ".params.rdata", sep = "")
            if (!dir.exists(label)) 
                dir.create(label)
            save(params, file = file.path(label, file))
            results$timing$end.time <- Sys.time()
            if (!is.null(num.secs)) {
                elapsed <- results$timing$end.time - results$timing$start.time
                units(elapsed) <- "secs"
                if (elapsed > num.secs) 
                  break
            }
        }
        elapsed <- results$timing$end.time - results$timing$start.time
        units(elapsed) <- "secs"
        completion <- num.iter * elapsed/i
        if (completion > 60) {
            units(completion) <- "mins"
            if (completion > 60) {
                units(completion) <- "hours"
                if (completion > 24) {
                  units(completion) <- "days"
                  if (completion > 7) {
                    units(completion) <- "weeks"
                  }
                }
            }
        }
        results$timing$completion.time <- completion
        results$timing$pct.complete <- round(100 * i/num.iter, 
            1)
        results$params <- params
    }, finally = setwd(wd))
    results
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
