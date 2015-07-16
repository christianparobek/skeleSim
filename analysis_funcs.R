#####################################################################


function(params){

  # saving global variables
  curr_scn<-params@current.scenario
  curr_rep<-params@current.replicate
  num_loci<-params@scenarios[[curr_scn]]@num.loci
  num_reps<-params@num.reps
  num_pops<-params@scenarios[[curr_scn]]@num.pops

  # If analysis results is empty, the first analysis done creates the list to hold the data
  if(is.null(params@analysis.results)){

    scenario.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)


    # creates a list of length scenarios for each requested groups of analyses
    for(group in names(which(params@analyses.requested)))
      scenario.results[[group]] <- vector('list', length(params@scenarios))

  } else {

    for(group in names(which(params@analyses.requested))){


      # group will cycle through selected among Global, Locus, and Pairwise
      # in each iteration of the for loop, group will have only one value

      if(group == "Global"){

        # TO DO check the data type and do conversions for what is needed ####@$@#$@#$#@

        # For multidna class objects we convert to a gtypes and use strataG for analysis
        #initialize arrays
        if (class(params@analysis.result)=="multidna"){

          num_loci <- getNumLoci(params@analysis.result[[2]])

          # Convert the list of DNAbin objects to gtypes
          genes <- params@analysis.result[[2]] #the multidna object
          names(genes@dna) <- paste("gene", 1:length(genes@dna))
          id <- genes@labels
          df <- data.frame(id = id, strata = params@analysis.result[[1]])
          gene.labels <- matrix(id, nrow = length(id), ncol = num_loci)
          colnames(gene.labels) <- paste("gene", 1:ncol(gene.labels), sep = "_")
          df <- cbind(df, gene.labels)
          results_gtype <- df2gtypes(df, 1, sequences = genes)


          #put overall analysis in first row using overall_stats()
          # params@current.replicate tells us how deep to put each new run in which list (@current.scenario)
          #run by locus analysis
          mat <- t(sapply(locNames(results_gtype), function (l){
            gtypes_this_loc<-subset(results_gtype, loci=l)
            overall_stats(gtypes_this_loc)
          }))
          analyses <- colnames(mat)
          num_analyses <- length(analyses)


          if(is.null(scenario.results[[group]][[curr_scn]])){
            scenario.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci,num_analyses,num_reps),
                                                           dimnames = list(1:num_loci,analyses,1:num_reps))
          }

          # We are printing by gene, not overall gene analysis. This differs from the genind code below.
          scenario.results[[group]][[curr_scn]][,,curr_rep] <- mat

          params@analysis.results <- scenario.results


        } else if(class(params@analysis.result)=="genind"){

          results_genind<-params@rep.result
          #convert genind to gtypes
          results_gtype<-genind2gtypes(results_genind)

          #   analyses <- names(overall_stats(results_gtype))
          #    num_analyses <- length(analyses)



          #put overall analysis in first row using overall_stats()
          # params@current.replicate tells us how deep to put each new run in which list (@current.scenario)
          #run by locus analysis
          mat <- t(sapply(locNames(results_gtype), function (l){
            gtypes_this_loc<-subset(results_gtype, loci=l)
            overall_stats(gtypes_this_loc)
          }))
          analyses <- colnames(mat)
          num_analyses <- length(analyses)

          # The first row will hold summary statistics over all loci regardless of population structure.
          # The remaining rows will hold summary statistics per locus
          if(is.null(scenario.results[[group]][[curr_scn]])){
            scenario.results[[group]][[curr_scn]] <- array(0, dim=c(1+num_loci,num_analyses,num_reps),
                                                           dimnames = list(c("Across_loci",1:num_loci),analyses,1:num_reps))
          }


          # combining overall statistics and locus by locus matrix
          scenario.results[[group]][[curr_scn]][,,curr_rep] <-   rbind(overall_stats(results_gtype),mat)

          params@analysis.results <- scenario.results

        }


        ################################### LOCUS ######################################
      } else if(x == "Locus"){

        if(class(params@analysis.result)=="genind"){

          results_genind<-params@rep.result

          #Hardy Weiberg test per population and overall
          pops_as_list<-seppop(results_genind)
          hw_results<-sapply(pops_as_list, function(p) hw.test(p)[,2:3], simplify = FALSE)
          hw_results.all <- rbind(hw.test(results_genind)[,2:3],do.call(rbind, hw_results)) # warnings for unknown reason
          colnames(hw_results.all)<-c("HWE.df","HWE.pval")


          #mratio on genind object, function needs genetic data as a genind, the population names, the locus names

          ##### Junk in here to get rid of!!!!
          pop.locus.df <- as.matrix(expand.grid(pop = levels(results_genind@pop), locus = levels(results_genind@loc.fac)))

          #  pop.locus.df <- as.matrix(expand.grid(pop = levels(nancycats@pop), locus = levels(nancycats@loc.fac)))

          mratio <- apply(pop.locus.df,1, function(x){
            calc.mratio(nancycats, x[1], x[2])
          })
          gtypes.cats <- genind2gtypes(nancycats)
          alleleFreqs(gtypes.cats, by.strata = TRUE)


          mratio <- apply(pop.locus.df,1, function(x){
            calc.mratio(results_genind, x[1], x[2])
          })
          
          #function to calculate Garza Williamson M ratio (bottleneck) statistic, per pop, per loc
          calc.gw <- function(gen.data.gtype) {
            pop.locus.df <- as.matrix(expand.grid(pop = 1:length(unique(gen.data.gtype@strata)), 
                                                  locus = 1:length(gen.data.gtype@loci)))
            
            all_freq_list<-alleleFreqs(gen.data.gtype, by.strata=T)
            mrat_pop_and_loc<-apply(pop.locus.df,1, function(x) {
              calc.mratio.2(all_freq_list[[x[2]]][,,x[1]][,1])              
            })
            
            calc.mratio.2<- function(all.freq){
              #find the smallest allele present
              if (any(all.freq >0)){ 
                min.all <- min(which(all.freq>0))
                #find the largest allele present
                max.all <- max(which(all.freq>0))
                #calculate the number of alleles present
                sum.all.present <- sum(all.freq>0)
                sum.all.present/(max.all-min.all+1)
              }
              else return(NA)
            }
            
            all_freq_list<-alleleFreqs(gen.data.gtype, by.strata=F)
            mrat_loc<-sapply(1:length(all_freq_list), function(x){
              calc.mratio.2(as.vector(all_freq_list[[x]][,1]))
            })
            return(c(mrat_loc,mrat_pop_and_loc))
            
          }


          #convert genind to gtypes for remaining analyses
          results_gtype<-genind2gtypes(results_genind)

          #by locus
          smryLoci <- cbind(Locus = 1:num_loci,Pop = NA,summarizeLoci(results_gtype))

          #by population
          locus_pop <- as.matrix(expand.grid(1:num_loci,1:num_pops))  #cycles through the loci for each population
          smryPop <- cbind(locus_pop,do.call(rbind,summarizeLoci(msats, by.strata = TRUE)))

          #Number of private alleles by locus
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


          ############### #mratio.p.val
          

         locus.final <- cbind(rbind(smryLoci,smryPop),hw_results.all)
         analysis_names <- colnames(locus.final)

          # Create the data array first time through
          if(is.null(scenario.results[[group]][[curr_scn]])){
            scenario.results[[group]][[curr_scn]] <- array(0, dim=c(num_loci*(num_pops+1),length(analysis_names),num_reps),
                                                           dimnames = list(1:(num_loci*(num_pops+1)),analysis_names,1:num_reps))

          }

          scenario.results[[group]][[curr_scn]][,,curr_rep] <-  locus.final

        }


        ############## ONCE ERIC fixes multidna stuffs, DO THIS!!!!  ###################
        if(class(params@analysis.result)=="multidna"){

        }

        # check location of this!
        params@analysis.results <- scenario.results


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