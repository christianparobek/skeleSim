analyses.check <- function(analyses.requested){
  #If no analyses are requested, default to all requested
  if(is.null(analyses.requested)){
    analyses.requested <- c(TRUE,TRUE,TRUE)
    names(analyses.requested) <- c("Global","Locus","Pairwise")
    analyses.requested
  } else {
    if(TRUE %in% analyses.requested){
      #if analyses exist but logicals are not named
      if(is.null(names(analyses.requested))){
        names(analyses.requested) <- c("Global","Locus","Pairwise")
        analyses.requested
      }
      } else {
        # no analyses.requested, default to all requested
        if(is.null(names(analyses.requested))){
          analyses.requested <- c(TRUE,TRUE,TRUE)
          names(analyses.requested) <- c("Global","Locus","Pairwise")
          analyses.requested
        } else {
          #no analyses.requested but they are named
          analyses.requested <- c(TRUE,TRUE,TRUE)
          analyses.requested
        }
      }
    }
  }
