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
  print(params@other.checks)
  #here we call the scenario checks (simulator specific and general)
  prv_chk <- params@sim.scen.checks  #store what is was in check slot
  #then calculate new checks
  ths_chk <- rbind(params@sim.check.func(params), gen.scenario.check(params))
  print(prv_chk);  print(ths_chk)
  #if what was there is null, replace with new checks
  if (is.null(prv_chk)) params@sim.scen.checks <- ths_chk
  #else, check which lines are there and replace info
  else {
    for (i in rownames(ths_chk)) {
      if (i %in% rownames(prv_chk)) prv_chk[i,] <- ths_chk[i,] #if it is there, replace it
      else {
        rbind(prv_chk,ths_chk[i,])  #if not, bind it
        rownames(prv_chk)[nrow(prv_chk)]<-i

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
    at.least.1.rep = params@num.reps > 0
  )
  print(results.check)
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

  results.check <- sapply(params@scenarios, function(sc) {
    #This will make a vector of TRUE/ FALSE
    c(nsizes.eq.npops = length(sc@pop.size) == sc@num.pops,
      nsamps.eq.npops = length(sc@sample.size) == sc@num.pops,
      is.mig.square = sapply(sc@migration, function(mig) {
        nrow(mig) == ncol(mig) & nrow(mig) == sc@num.pops
      }),
      at.lst.1.pop = sc@num.pops >= 1,
      at.lst.1.loc = all(sc@num.loci >= 1),
      mut.rate.ok = all((sc@mut.rate>=0)&(sc@mut.rate<1)),
      at.lst.1.samp = min(sc@sample.size)>0,
      mig.bet.0.1 = sapply(sc@migration, function(mig) {
        all((mig>=0)&(mig<=1))
      }),
      mig.diag.eq.0 = sapply(sc@migration, function(mig) {
        all(diag(mig)==0)
      })
    )
  })
  return(results.check)
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