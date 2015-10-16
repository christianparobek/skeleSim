#####################################################################


# TO DO: remove the scenario.results and assign the framework for repetitions to params@analysis.results
# analyse param@sample which is the

# To Do: not all the assignments to params@analysis.results are to the correct spot in the list, need to check

# To Do: once the if/else statements are changed to just if statements should remove an extra list creation
# for params@analysis.results... I think there are extra.

function(params){

  # saving global variables
  curr_scn<-params@current.scenario
  curr_rep<-params@current.replicate
  num_loci<-params@scenarios[[curr_scn]]@num.loci
  num_reps<-params@num.reps
  num_pops<-params@scenarios[[curr_scn]]@num.pops

  #repsample data - create a gtype then just switch to genind when need to

  # If analysis results is empty, the first analysis done creates the list to hold the data
  if(is.null(params@analysis.results)){

    params@analysis.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)

    # creates a list of length scenarios for each requested groups of analyses
    for(group in names(which(params@analyses.requested)))
      params@analysis.results[[group]] <- vector('list', length(params@scenarios))

# TO DO: get rid of this else.
} else {

    for(group in names(which(params@analyses.requested))){

      # group will cycle through selected among Global, Locus, and Pairwise
      # in each iteration of the for loop, group will have only one value

      if(group == "Global"){

        # TO DO check the data type and do conversions for what is needed
        # For multidna class objects we convert to a gtypes and use strataG for analysis

        # if genes > 1 do different formatting

        #initialize arrays
        if (class(params@analysis.results)=="multidna"){

          num_loci <- nLoc(params@rep.sample[[2]])

          # Convert the list of DNAbin objects to gtypes
          genes <- params@analysis.result[[2]] #the multidna object
          names(genes@dna) <- paste("gene", 1:length(genes@dna))
          id <- genes@labels
          df <- data.frame(id = id, strata = params@analysis.result[[1]])
          gene.labels <- matrix(id, nrow = length(id), ncol = num_loci)
          colnames(gene.labels) <- paste("gene", 1:ncol(gene.labels), sep = "_")
          df <- cbind(df, gene.labels)
          results_gtype <- df2gtypes(df, 1)


          results.matrix <- t(sapply(locNames(results_gtype), function (l){
            gtypes_this_loc<-results_gtype[,l,]
            overall_stats(gtypes_this_loc)
          }))
          analyses <- colnames(results.matrix)
          num_analyses <- length(analyses)


          if(is.null(params@analysis.results[[group]][[curr_scn]])){
            params@analysis.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci,num_analyses,num_reps),
                                                           dimnames = list(1:num_loci,analyses,1:num_reps))
          }

          #this shouldn't happen here it should be at close of function- this is what is returned
          # We are printing by gene, not overall gene analysis. This differs from the genind code below.
          params@analysis.results[[group]][[curr_scn]][,,curr_rep] <- results.matrix




        } else if(class(params@analysis.result)=="genind"){

          #Global

          results_genind<-params@rep.result
          #convert genind to gtypes
          results_gtype<-genind2gtypes(results_genind)

          #   analyses <- names(overall_stats(results_gtype))
          #    num_analyses <- length(analyses)

          #put overall analysis in first row using overall_stats()
          # params@current.replicate tells us how deep to put each new run in which list (@current.scenario)
          #run by locus analysis
          mat <- t(sapply(locNames(results_gtype), function (l){
            gtypes_this_loc<-results_gtype[,l,]
            overall_stats(gtypes_this_loc)
          }))
          analyses <- colnames(mat)
          num_analyses <- length(analyses)

          # The first row will hold summary statistics over all loci regardless of population structure.
          # The remaining rows will hold summary statistics per locus
          if(is.null(params@analysis.results[[group]][[curr_scn]])){
            params@analysis.results[[group]][[curr_scn]] <- array(0, dim=c(1+num_loci,num_analyses,num_reps),
                                                           dimnames = list(c("Across_loci",1:num_loci),analyses,1:num_reps))
          }

          #this shouldn't happen here it should be at close of function- this is what is returned
          # combining overall statistics and locus by locus matrix
          params@analysis.results[[group]][[curr_scn]][,,curr_rep] <-   rbind(overall_stats(results_gtype),mat)



        }


        ################################### LOCUS ######################################
      } else if(x == "Locus"){

        if(class(params@analysis.result)=="genind"){

          results_genind<-params@rep.result

          #Hardy Weiberg test per population and overall (this comes first because it needs genind)
          pops_as_list<-seppop(results_genind)
          hw_results<-sapply(pops_as_list, function(p) hw.test(p)[,2:3], simplify = FALSE)
          hw_results.all <- rbind(hw.test(results_genind)[,2:3],do.call(rbind, hw_results)) # warnings for unknown reason
          colnames(hw_results.all)<-c("HWE.df","HWE.pval")

          #convert genind to gtypes for remaining analyses
          results_gtype<-genind2gtypes(results_genind)

          #mratio on gtypes object, function needs genetic data as a gtype
          mrat_results_all<-calc.mratio(results_gtype)

          #by locus, all the other stats (num alleles etc) pulled from summarizeLoci
          smryLoci <- cbind(Locus = 1:num_loci,Pop = NA,summarizeLoci(results_gtype))

          #by population
          locus_pop <- as.matrix(expand.grid(1:num_loci,1:num_pops))  #cycles through the loci for each population
          smryPop <- cbind(locus_pop,do.call(rbind,summarizeLoci(results_gtype, by.strata = TRUE)))

          #Number of private alleles by locus- this should be moved to the skelesimFuncs
          alleleFreqs <- alleleFreqs(results_gtype, by.strata = TRUE)
          by.loc <- sapply(alleleFreqs, function(loc) {
            mat <- loc[, "freq", ]
            rowSums(apply(mat, 1, function(r) {
              result <- rep(FALSE, length(r))
              if(sum(r > 0) == 1) result[r > 0] <- TRUE
              result
            }))
          })
          rownames(by.loc) <- strataNames(results_gtype)
          perLocus <- colSums(by.loc) #this has the number of alleles that are private per locus

          #the rows will be have the private alleles for each population by locus
          smryPop <- cbind(smryPop, as.vector(t(by.loc)))
          smryLoci <- cbind(smryLoci, num.priv.all = perLocus)


            ################# #Convert from genind to loci (package pegas))
          #Allan, Sean and I decided that a workable alternative to Fis estimation below is to use:   Fis = 1-(Ho/He) since we already have obs.He and exp.He is the results
          results_loci<-genind2loci(results_genind)
          #for loci
          FSTloci<-Fst(results_loci)
          Fisloci <- FSTloci[ , c("Fis")]
          Fitloci <- FSTloci[ , c("Fit")]
          Fstloci <- FSTloci[ , c("Fst")]
          #for pops
          FSTpop<-Fst(results_loci, pop = results_loci$population) #column with pop information#)
                        Fispop <- FSTpop[ , c("Fis")]
                      Fitpop <- FSTpop[ , c("Fit")]
                      Fstpop <- FSTpop[ , c("Fst")]


          #all the analyses get bound here
         locus.final <- cbind(rbind(smryLoci,smryPop),hw_results.all,mrat_results_all)
         analysis_names <- colnames(locus.final)

          # Create the data array first time through
          if(is.null(params@analysis.results[[group]][[curr_scn]])){
            params@analysis.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci*(num_pops+1),length(analysis_names),num_reps),
                                                           dimnames = list(1:(num_loci*(num_pops+1)),analysis_names,1:num_reps))

          }

          params@analysis.results[[group]][[curr_scn]][,,curr_rep] <-  locus.final

        }


        ############## ONCE ERIC fixes multidna stuffs, DO THIS!!!!  ###################
        if(class(params@analysis.result)=="multidna"){

          overallTest




        }


###########################  Pairwise  ###########################
      } else if(x == "Pairwise"){

        ##### no difference between pws and genotype data pws
        #Pairwise Chi2, D, F..., G...

        #dA, mean.pct.between


        ###Should be able to remove these two lines for the next block to deal with any possible multiDNA objects
        #pws <- pairwiseTest(results_gtype, nrep = 5,stat.list = list(statGst), quietly = TRUE)
        #pws[1]$result[,-c(2:5)] # assuming the 2nd through 5th column remain strata.1, strata.2, n.1, and n.2

        #Pairwise Chi2, Fst, PHIst
        pws.multi <- sapply(locNames(results_gtype), function(g){
          gene_gtype <- results_gtype[,g,]
          pairwiseTest(gene_gtype, nrep =5, keep.null=TRUE)
        })
        pws.mulit[1]$result[,-c(2:5)] #removing strata.1, strata.2, n.1, n.2

        #Genotype data
        pws <- pairwiseTest(results_gtype, nrep = 5, stat.list = list(statGst, quietly = TRUE))
        pws.out <- pws$result[-c(2:5)]
        # rownames(pws.out) <- as.matrix(pws[[1]][1])  # could turn into row names at the end.. but they'll be repeated
        sA <- sharedAlleles(results_gtype)[,-c(1:2)]
        nsharedAlleles <- paste("sharedAlleles", names(sA), sep = ".")
        names(sA) <- nsharedAlleles
        pws.out <- cbind(pws.out, sA)
        params@analysis.results[[curr_scn]] <- pws.outS
      }

      #Data.frame of summary data into simulation replicate

      # Enter array of data into the current scenario
      params@analysis.results[[curr_scn]] <- data.cube


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
  rep.sample.gtypes <- df2gtypes(df, 1)

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