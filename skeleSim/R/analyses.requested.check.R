#' @title Ensure Analyses Request Vector is Well-formed
#' @description take a named vector and make sure that it contains the names Global, Locus, and Pairwise
#'
#' @param analyses.requested A named logical vector with elements named Global, Locus, and Pairwise
#'
#' @export
analyses.check <- function(analyses.requested)
{
    ret <- NULL
    if (sum(!names(analyses.requested)%in%c("Global","Locus","Pairwise"))>0)
    {
        stop("there are analyses requested that we dont recognize")
    }
    if (length(names(analyses.requested))==0) {analyses.requested=NULL} #if the vector of inputs is not named, assume all are requested
    if (is.null(analyses.requested))
        {
            ret <- c("Global"=T,"Pairwise"=T,"Locus"=T)
        } else {

            add <- c("Global","Locus","Pairwise")[!(c("Global","Locus","Pairwise")%in%names(analyses.requested))]

            if (length(add)>0)
            {
                tlen <- length(analyses.requested)
                analyses.requested <- c(analyses.requested,rep(F,length(add)))
                names(analyses.requested)[names(analyses.requested)==""] <- add
            }
            ret <- analyses.requested
        }
    ret
}

