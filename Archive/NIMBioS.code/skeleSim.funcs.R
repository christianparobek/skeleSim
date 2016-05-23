currentLabel <- function(params) {
  label <- paste(
    params@title,
    params@current.scenario,
    params@current.replicate,
    sep = "."
  )
  gsub("[[:punct:]]", ".", label)
}


currentScenario <- function(params) {
  params@scenarios[[params@current.scenario]]
}


overall.check <- function(params) {
  params@other.checks <- non.scenario.check(params)
  print(params@other.checks)
  #here we call the scenario checks (simulator specific and general)
  prv_chk<-params@sim.scen.checks  #store what is was in check slot
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


non.scenario.check<- function(params) {
  results.check <- c(
    title.not.null = !is.null(params@title),
    #check that number of reps is greater than 0
    at.least.1.rep = params@num.reps > 0
  )
  print(results.check)
  return(results.check)
}

gen.scenario.check <-function(params) {
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
      at.lst.1.loc = sc@num.loci >= 1,
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


runSim <- function(params, num.secs = NULL) {
  # check parameters
  params<-overall.check(params)
  if(!all(params@other.checks, params@sim.scen.checks)) {
    filewr <- "error.log"
    write.csv(params@sim.scen.checks, file=filewr)
    print(params@other.checks)
    write("\n",file=filewr,append=T)
    write.table(params@other.checks, file=filewr,append=T)

    stop("parameters do not pass checks; see error log for details")
  }
print("params checked")
  
  # Check/setup folder structure
  if(file.exists(params@wd)) {
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
    # loop through replicates for scenarios
    num.iter <- nrow(params@scenario.reps)
    results$timing$start.time <- Sys.time()
    for(i in 1:num.iter) {
      if(quit) break # leave replicate for loop if user has decided to quit
      params@current.scenario <- params@scenario.reps[i, "scenario"]
      params@current.replicate <- params@scenario.reps[i, "replicate"]
      # run one replicate of simulator
      params <- params@sim.func(params)
      # analyzes params@rep.sample and loads results into params@rep.result
      params <- params@rep.analysis.func(params)
# -->> REMOVE FOR RELEASE: SAVING params OBJECT FOR TESTING <<--
      label <- currentLabel(params)
      file <- paste(label, ".params.rdata", sep = "")
      if(!dir.exists(label)) dir.create(label)
      save(params, file = file.path(label, file))
      # check timing
      results$timing$end.time <- Sys.time()
      if(!is.null(num.secs)) {
        elapsed <- results$timing$end.time - results$timing$start.time
        units(elapsed) <- "secs"
        if(elapsed > num.secs) break
      }
    }

    elapsed <- results$timing$end.time - results$timing$start.time
    units(elapsed) <- "secs"
    completion <- num.iter * elapsed / i
    if(completion > 60) {
      units(completion) <- "mins"
      if(completion > 60) {
        units(completion) <- "hours"
        if(completion > 24) {
          units(completion) <- "days"
          if(completion > 7) {
            units(completion) <- "weeks"
          }
        }
      }
    }
    results$timing$completion.time <- completion
    results$timing$pct.complete <- round(100 * i / num.iter, 1)
    results$params <- params
  }, finally = setwd(wd))

  results
}


summ.stats.table<-function(params){

  results.datafr<-params@analysis.results
  num.sc <- length(params@scenarios)
  names.stats<-colnames(params@analysis.results)
  colnames(results.datafr)<-c(names.stats)
  num.stats<-length(names.stats)-1

  #table of means and SD
  table.means<-data.frame(matrix(NA,nrow=num.sc,ncol=num.stats))
  table.sd<-data.frame(matrix(NA,nrow=num.sc,ncol=num.stats))
  colnames(table.means)<-names.stats[-1]; colnames(table.sd)<-names.stats[-1]
  for (i in 1:num.stats)
    table.means[,i]<-tapply(results.datafr[,1+i],results.datafr[,1],mean)
  for (i in 1:num.stats)
    table.sd[,i]<-tapply(results.datafr[,1+i],results.datafr[,1],sd)
  results.list<-list()
  results.list$means<-table.means
  results.list$sd<-table.sd
  params@summary.results<-results.list

  return(params)
}


plot.all.stats<-function(params){
  stopifnot(require("reshape2"))
  stopifnot(require("ggplot2"))
  stopifnot(require("gridExtra"))

  results.datafr<-as.data.frame(params@analysis.results)
  num.sc <- length(params@scenarios)
  names.stats<-colnames(params@analysis.results)
  colnames(results.datafr)<-c(names.stats)
  num.stats<-length(names.stats)-1

  #plotting option 1
  results.melted<-melt(results.datafr,id="scenario",measure.vars=c(names.stats[-1]))
  ggplot(results.melted, aes(value)) +
      geom_density() + facet_wrap(scenario~variable,
      ncol = num.stats, scales = "free")

  #histogram plot
  ggplot(results.melted, aes(value)) + geom_histogram() +
  facet_grid(variable~scenario, scales="free")

 }

#run overall analysis
overall_stats<- function(results_gtype){
  ovl <- overallTest(results_gtype, nrep = 5, quietly = TRUE)
  ovl.result <- ovl$result[complete.cases(ovl$result[,1]),]

  pnam <- c()
  for(i in 1:nrow(ovl.result))
    pnam <- c(pnam,rownames(ovl.result)[i],paste(rownames(ovl.result)[i],"pval", sep = ""))

  global.wide <- as.vector(t(ovl.result))
  names(global.wide) <- pnam
  global.wide
}
