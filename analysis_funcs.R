#####################################################################
####      ROW NAMES
# for pairwise populations, take number of populations = 'params@num.pops'

# params@rep.sample
# params@num.pops


### create containers with different sizes by scenario
# as many cubes as scenario, the size is the replicate
# rows will be number of loci and this will change/vary by scenario
# in order to look at number of loci needed (power analysis)
# loop for x in which(names(params@analysis.requested))

# params@rep.result will be genind, after can overwrite @rep.result with gtypes, can assign NULL to @rep.results
# simulation -> genind - do all calcs -> gtypes conversion - do all calculations (overwrite, only one object)

####      ROW NAMES
# for pairwise populations, take number of populations = 'npop'
npp <- combn(1:npop, 2)
names <- c("overall", apply(npp, 2, function(x) paste(x, collapse = "v")))

# for pairwise loci, take number of loci = nloc
npl <- combn(1:nloc, 2)
names <- c("overall", apply(npl, 2, function(x) paste(x, collapse = "v")))

# for simple populations
names <-  c("overall", 1:npop)

# for loci
names <- c("overall", 1:nloc)
#####################################################################


function(params){

  curr_scn<-params@current.scenario
  curr_rep<-params@current.replicate
  num_loci<-params@scenarios[[curr_scn]]@num.loci
  num_reps<-params@num.reps
  num_pops<-params@scenarios[[curr_scn]]@num.pops


if(is.null(params@analysis.results)){

  scenario.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)



  for(group in names(which(params@analyses.requested)))
    scenario.results[[x]] <- vector('list', length(params@scenarios))

} else {

  for(group in names(which(params@analyses.requested))){


    # group will cycle through selected among Global, Popualtion, Locus, and Pairwise
    # in each iteration of the for loop, group will have only one value!

if(group == "Global"){

  # TO DO check the data type and do conversions for what is needed

  #initialize arrays
  if (class(params@analysis.result)=="multidna"){

    analyses <- names(overall_stats(results_gtype))
    num_analyses <- length(analyses)

    genes <- params@analysis.result[[2]]
    names(genes@dna) <- paste("gene", 1:length(genes@dna))
    id <- genes@labels
    df <- data.frame(id = id, strata = params@analysis.result[[1]], hap = id)
    results_gtype<-df2gtypes(df, 1, sequences = genes)


      scenario.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci,num_analyses,num_reps),
                                                     dimnames = list(1:num_loci,analyses,1:num_reps))


################################### once Eric fixes naming the loci, start here for multidna objects global for each gene ###########



    } else if(class(params@analysis.result)=="genind"){

    results_genind<-params@rep.result
    #convert genind to gtypes
    results_gtype<-genind2gtypes(results_genind)

    analyses <- names(overall_stats(results_gtype))
    num_analyses <- length(analyses)


    # The first row will hold summary statistics over all loci regardless of population structure.
    # The remaining rows will hold summary statistics per locus
    scenario.results[[group]][[curr_scn]] <- array(0, dim=c(1+num_loci,num_analyses,num_reps),
                                                  dimnames = list(c("Across_loci",1:num_loci),analyses,1:num_reps))

  #put overall analysis in first row using overall_stats()
  # params@current.replicate tells us how deep to put each new run in which list (@current.scenario)
  #run by locus analysis
  mat <- t(sapply(locNames(results_gtype), function (l){
    gtypes_this_loc<-subset(results_gtype, loci=l)
    overall_stats(gtypes_this_loc)
  }))

  # combining overall statistics and locus by locus matrix
  scenario.results[[group]][[curr_scn]][,,curr_rep] <-   rbind(overall_stats(results_gtype),mat)

  params@analysis.results <- scenario.results

}


################################### LOCUS ######################################
} else if(x == "Locus"){

  else if(class(params@analysis.result)=="genind"){

    # genind to gtypes - split -> genind of split and then do

    results_genind<-params@rep.result

    seppop(results_genind)

    #Hardy Weiberg test
    pops_as_list<-seppop(results_genind)
    hw_results<-sapply(pops_as_list, function(p) hw.test(p)[,2:3], simplify = FALSE)
    hw_results.all <- rbind(hw.test(results_genind)[,2:3],do.call(rbind, hw_results))
    analysis_names2 <- colnames(hw_results.all)


    #convert genind to gtypes
    results_gtype<-genind2gtypes(results_genind)

 # analysis_names <- colnames(summarizeLoci(results_gtype))
  locus_pop <- as.matrix(expand.grid(1:num_loci,1:num_pops))  #cycles through the loci for each population

 #by loci
  smryLoci <- cbind(Locus = 1:num_loci,Pop = NA,summarizeLoci(msats))
  analysis_names1 <- colnames(smryLoci)

 #by population
  smryPop <- cbind(locus_pop,do.call(rbind,summarizeLoci(msats, by.strata = TRUE)))

 if(is.null(scenario.results[[group]][[curr_scn]])){
   scenario.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci*(num_pops+1),length(analysis_names)+4,num_reps),
                                                dimnames = list(1:(num_loci*(num_pops+1)),c(analysis_names1,analysis_names2),1:num_reps))

 }

  scenario.results[[group]][[curr_scn]][,,curr_rep] <-  cbind(rbind(smryLoci,smryPop),hw_results.all)

  }


} else if(x == "Pairwise"){
  #Genotype data
  pws <- pairwiseTest(params@rep.result, nrep = 5, stat.list = list(statGst, quietly = TRUE))
  pws.out <- pws$result[-c(1:5)]
  rownames(pws.out) <- as.matrix(pws[[1]][1])
  sA <- sharedAlleles(params@rep.result)[,-c(1:2)]
  nsharedAlleles <- paste("sharedAlleles", names(sA), sep = ".")
  names(sA) <- nsharedAlleles
  pws.out <- cbind(pws.out, sA)
    scenario.results[[params@current.scenario]] <- pws.out
}

  #Data.frame of summary data into simulation replicate

  # Enter array of data into the current scenario
  scenario.results[[params@current.scenario]] <- data.cube


  }

}

params

}


##### Will need to be in loop
  # params@num.pops <- the number of rows, pops or pairwise or loci plus overall
  # nrep <- the number of replicates
  # dimension names will come from the choice of pairwise, by
  #   by population, by loci, or whatever


  sim.data.analysis <- array(0, dim = c(length(names), length(analyses), nrep),
                             dimnames = list(names, analyses, 1:nrep))




  ### loading function
  # This gets called and knows which row (so z spot) matchs the current
  #   scenario and replicate
  ########### multiDNA to gtypes
  if(class(params@rep.sample) == "multidna"){
    genes <- rep.sample[[2]]
    names(genes@dna) <- paste("gene", 1:length(genes@dna))
    id <- genes@labels
    df <- data.frame(id = id, strata = rep.sample[[1]], hap = id)
    rep.sample.gtypes <- df2gtypes(df, 1, sequences = genes)

  }

  if(class(sim.out) == "DNAbin"

     if(pairwisepop = TRUE){

       npp <- combn(1:params@num.pops, 2)
       names <- c("overall", apply(npp, 2, function(x) paste(x, collapse = "v")))
     } else if (pairwiselocus = TRUE){
       # for pairwise loci, take number of loci = nloc
       npl <- combn(1:nloc, 2)
       names <- c("overall", apply(npl, 2, function(x) paste(x, collapse = "v")))
     } else if(bypop = TRUE){
       # for simple populations
       names <-  c("overall", 1:params@num.pops)
     } else {
       # for loci
       names <- c("overall", 1:nloc)
     }
     #####################################################################

     # testing ones
     analyses <- c("allel","Freq","prop")
     nrep <- 5


     #####################################################################

     # "all.data" needs to be merged into one matrix per replicate
     ### Global, Genotypes
     # Chi2, D, Fst, F'st, Gst, G'st, G''st, P-vals for each
     ovl <- overallTest(simdata, nrep = 5, stat.list = statList("chi2"), quietly = TRUE)
     global <- t(ovl$result)
     pnam <- c()
     for(i in 1:length(colnames(global))){
       pnam <- c(pnam,paste(colnames(global)[i],"pval", sep = ""))
     }
     global.wide <- c(global[1,],global[2,])
     names(global.wide) <- c(colnames(global),pnam)


     ############################################################################

     for(i in 1:params@num.reps){
       sim.data.analysis[,,i] <- matrix(all.data)
     }

     ###########################################################################