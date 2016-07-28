#' @name skeleSim.checks
#' @title Check all simulation parameters
#' @description Check all simulation parameters
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @export
#'
overall.check <- function(params) {
  params@other.checks <- non.scenario.check(params)
  #here we call the scenario checks (simulator specific and general)
  prv_chk <- params@sim.scen.checks  #store what is was in check slot
  #then calculate new checks
  ths_chk <- rbind(params@sim.check.func(params), gen.scenario.check(params))

  #if the number of columns of the two checks are different, then set prv_chk to null
  #this would happen if a scenario is added to an existing object
  if (!is.null(prv_chk))
      if (dim(prv_chk)[2]!=dim(ths_chk)[2])
          prv_chk=NULL

  
  if(is.null(prv_chk)) { #if what was there is null, replace with new checks
    params@sim.scen.checks <- ths_chk
  } else { # else, check which lines are there and replace info
    for(x in rownames(ths_chk)) {
      if(x %in% rownames(prv_chk)) { #if it is there, replace it
        prv_chk[x, ] <- ths_chk[x, ]
      } else { # if not, bind it
        prv_chk <- rbind(prv_chk, ths_chk[x, ])
        rownames(prv_chk)[nrow(prv_chk)] <- x
      }
    }
    params@sim.scen.checks <- prv_chk
  }

  ##############TO DO write to a file error log#################

  #output result based on both sets of checks
  return(params)
}


#' @rdname skeleSim.checks
#' @export
#'
non.scenario.check <- function(params) {
  results.check <- c(
    title.not.null = !is.null(params@title),
    #check that number of reps is greater than 0
    at.least.1.rep = params@num.sim.reps > 0
  )
  return(results.check)
}


#' @rdname skeleSim.checks
#' @export
#'
gen.scenario.check <- function(params) {
  #check that number of populations is same as length of pop sizes
  #check that number of populations is same as length of sample sizes
  #check that migration matrix is square with sises equal to number pops
  #check that num.pops has to be 1 or greater
  #check that num.loci has to be 1 or greater
  #check that mut.rate between 0 (inclusive) and 1 (exclusive)
  #check that at least one sample sizes has to be greater than 0
  #check that migration matrix all between 0 and 1
  #check that all migration matrix diagonals are 0

  sapply(params@scenarios, function(sc) {
    mig <- sc@migration
    mig.fmt <- if(is.null(mig)) TRUE else {
      if(is.matrix(mig)) mig <- list(mig)
      if(is.list(mig)) all(sapply(mig, is.matrix)) else FALSE
    }
    is.mig.square <- if(!mig.fmt) FALSE else {
      if(is.null(mig)) TRUE else {
        all(sapply(mig, function(x) {
          nrow(x) == ncol(x) & nrow(x) == sc@num.pops
        }))
      }
    }
    mig.btwn.0.1 <- if(!mig.fmt) FALSE else {
      if(is.null(mig)) TRUE else all(sapply(mig, function(x) all(x >= 0 & x <= 1)))
    }
    mig.diag.eq.0 <- if(!mig.fmt) FALSE else {
      if(is.null(mig)) TRUE else all(sapply(mig, function(x) all(diag(x) == 0)))
    }

    c(nsizes.eq.npops = length(sc@pop.size) == sc@num.pops,
      nsamps.eq.npops = length(sc@sample.size) == sc@num.pops,
      at.lst.1.pop = sc@num.pops >= 1,
      at.lst.1.loc = all(sc@num.loci >= 1),
      at.lst.1.samp = min(sc@sample.size) > 0,
      mut.rate.ok = all(sc@mut.rate >= 0 & sc@mut.rate < 1),
      mig.fmt = mig.fmt,
      is.mig.square = is.mig.square,
      mig.btwn.0.1 = mig.btwn.0.1,
      mig.diag.eq.0 = mig.diag.eq.0
    )
  })
}


#' @rdname skeleSim.checks
#' @param analyses.requested A named logical vector with elements named Global, Locus, and Pairwise
#'
#' @export
#'
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
