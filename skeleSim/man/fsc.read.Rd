\name{fsc.read}
\alias{fsc.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fsc.read(file, sc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
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
function (file, sc) 
{
    stopifnot(require(adegenet) & require(apex))
    f <- readLines(file)
    start <- grep("SampleData=", f) + 1
    end <- which(f == "}") - 2
    pos <- cbind(start, end)
    pop.data <- do.call(rbind, lapply(1:nrow(pos), function(i) {
        f.line <- f[pos[i, 1]:pos[i, 2]]
        f.line <- gsub("[[:space:]]+", "--", f.line)
        data.mat <- do.call(rbind, strsplit(f.line, "--"))[, 
            -2]
        data.mat <- cbind(rep(paste("Sample", i), nrow(data.mat)), 
            data.mat)
    }))
    data.type <- f[grep("DataType=", f)]
    data.type <- gsub("\tDataType=", "", data.type)
    is.seq <- switch(data.type, DNA = T, MICROSAT = F, STANDARD = F)
    if (is.seq) {
        poly.pos.lines <- grep("polymorphic positions on chromosome", 
            f, value = T)
        num.poly <- as.numeric(sapply(strsplit(poly.pos.lines, 
            " "), function(x) x[2]))
        end <- cumsum(num.poly)
        start <- c(1, end[1:(length(end) - 1)] + 1)
        seq.len <- sc@sequence.length
        dna.seqs <- lapply(1:length(seq.len), function(i) {
            seq.i <- if (num.poly[i] == 0) {
                full.seq <- paste(rep("A", seq.len[i]), collapse = "")
                rep(full.seq, nrow(pop.data))
            }
            else {
                padding <- paste(rep("A", seq.len[i] - num.poly[i]), 
                  collapse = "")
                seq.i <- substr(pop.data[, 3], start[i], end[i])
                sapply(seq.i, function(x) paste(x, padding, sep = "", 
                  collapse = ""))
            }
            names(seq.i) <- pop.data[, 2]
            as.DNAbin(strsplit(tolower(seq.i), ""))
        })
        dna.seqs <- new("multidna", dna.seqs)
        list(strata = pop.data[, 1], dna.seqs = dna.seqs)
    }
    else {
        n.loc <- ncol(pop.data) - 2
        pop.data <- do.call(rbind, lapply(seq(1, nrow(pop.data), 
            2), function(i) {
            ind <- pop.data[c(i, i + 1), ]
            locus.data <- as.vector(ind[, -(1:2)])
            c(ind[1, 1], paste(ind[, 2], collapse = "/"), locus.data)
        }))
        locus.data <- pop.data[, -c(1:2)]
        collapsed.loci <- do.call(cbind, lapply(seq(2, ncol(locus.data), 
            by = 2), function(i) {
            a1 <- locus.data[, i - 1]
            a2 <- locus.data[, i]
            paste(a1, a2, sep = "/")
        }))
        colnames(collapsed.loci) <- paste("Locus", 1:ncol(collapsed.loci), 
            sep = ".")
        rownames(collapsed.loci) <- pop.data[, 2]
        df2genind(collapsed.loci, sep = "/", pop = pop.data[, 
            1], type = "codom")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
