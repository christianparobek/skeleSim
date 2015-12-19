function(params){

  stopifnot(require(strataG))
  stopifnot(require(pegas))

  # saving global variables
  curr_scn<-params@current.scenario
  curr_rep<-params@current.replicate
  num_loci<-params@scenarios[[curr_scn]]@num.loci
  num_reps<-params@num.reps
  num_pops<-params@scenarios[[curr_scn]]@num.pops

  #params@rep.sample is either a genind or a list of DNAbin objects
  if(inherits(params@rep.sample, "genind")){
    results_gtype<-genind2gtypes(params@rep.sample)
    } else if (inherits(params@rep.sample, "gtypes")){  #should -> "DNAbin"? and first do new("multidna", list(sequences)) and sequence2gtypes(new("multidna"...))
      results_gtype <- params@rep.sample
    } else {
      # Convert the list of DNAbin objects to gtypes
      genes <- params@rep.sample$dna.seqs #the multidna object
      names(genes@dna) <- paste("gene", 1:length(genes@dna))
      id <- genes@labels
      df <- data.frame(id = id, strata = params@rep.sample$strata)
      gene.labels <- matrix(id, nrow = length(id), ncol = num_loci)
      colnames(gene.labels) <- paste("gene", 1:ncol(gene.labels), sep = "_")
      df <- cbind(df, gene.labels)
      results_gtype <- df2gtypes(df, 1)
    }

  # If analysis results is empty, the first analysis done creates the list to hold the data
  if(is.null(params@analysis.results)){
    params@analysis.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)
  }

        ######################### Global ##########################

  # TO DO: Need an across all loci row if a genind object in addition to each loci Will need a different
  # analysis.results format

  if(params@analyses.requested["Global"]){

    #Don't need this, it's a global variable
    # num_loci <- nLoc(results_gtype)

    #overall_stats() is from skeleSim.funcs.R
    #run by locus analysis across all populations
    r.m <- lapply(locNames(results_gtype), function (l){
      gtypes_this_loc<-results_gtype[,l,]
      overall_stats(gtypes_this_loc)
      })
    results.matrix <- do.call(rbind, r.m)
    analyses <- colnames(results.matrix)
    num_analyses <- length(analyses)
    rownames(results.matrix) <- locNames(results_gtype)

    # We are printing by gene, not overall gene analysis. This differs from the genind code below.
    # The first row will hold summary statistics over all loci regardless of population structure.
    # The remaining rows will hold summary statistics per locus
    if(is.null(params@analysis.results[["Global"]][[curr_scn]])){
      params@analysis.results[["Global"]][[curr_scn]] <- array(0,dim=c(num_loci,
                                                                       num_analyses,
                                                                       num_reps),
                                                               dimnames = list(c(1:num_loci),
                                                                               analyses,
                                                                               1:num_reps))
    }

    params@analysis.results[["Global"]][[curr_scn]][,,curr_rep] <- results.matrix
  }
}


        ################################### LOCUS ######################################

  if(params@analyses.requested["Locus"]){

# genind
        if(inherits(params@rep.sample,"genind")){

          results_genind<-params@rep.sample

          #Hardy Weiberg test per population and overall (this comes first because it needs genind)
          pops_as_list<-seppop(results_genind)
          hw_results<-sapply(pops_as_list, function(p) hw.test(p)[,2:3], simplify = FALSE)
          hw_results.all <- rbind(hw.test(results_genind)[,2:3],do.call(rbind, hw_results))
          colnames(hw_results.all)<-c("HWE.df","HWE.pval")

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


          ################# Convert from genind to loci (package pegas))
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
          if(is.null(params@analysis.results[["Locus"]][[curr_scn]])){
            params@analysis.results[["Locus"]][[curr_scn]] <- array(0, dim=c(num_loci*(num_pops+1),
                                                                             length(analysis_names),
                                                                             num_reps),
                                                                    dimnames = list(1:(num_loci*(num_pops+1)),
                                                                                    analysis_names,
                                                                                    1:num_reps))
          }

          params@analysis.results[["Locus"]][[curr_scn]][,,curr_rep] <-  locus.final

        }

    # multiDNA
    # Per gene (ignoring population structure, and per gene by population)

      if(inherits(params@rep.sample,c("multidna","gtypes"))){

          # by gene
          r.m.gene <- lapply(locNames(results_gtype), function(l){
            nucleotideDiversity(results_gtype[,l,]@sequences)
          })
          results.list.names.gene <- Map(function(gene,names) paste(gene, names(names), sep="_"),
                                         locNames(results_gtype),
                                         r.m.gene)
          results.gene <- do.call(c,r.m.gene)
          names(results.gene) <- do.call(c,results.list.names.gene)

          # by gene per popualation "strata"
          r.m <- lapply(strataNames(results_gtype), function(s){
            lapply(locNames(results_gtype), function(l){
            nucleotideDiversity(results_gtype[,l,s]@sequences)
            })
          })

          r.m.bind <- do.call(c, lapply(r.m, function(x){
            do.call(c,x)
          }))

          results.list.names <- lapply(1:length(strataNames(results_gtype)), function(x){
            Map(function(gene,names) paste(gene, names(names), sep="_"),
                locNames(results_gtype),
                r.m[[x]])
          })

          rln.bind <- do.call(c, lapply(results.list.names, function(x){
            do.call(c,x)
             }))

          results <- r.m.bind
          #names(results) <- rln.bind

          nD <- rbind(results.gene, results)


          #fusFs
          ### Warning: Some sequences could not be unambiguously assigned to a haplotype
          #by gene
          fu.fs <- lapply(locNames(results_gtype), function(l){
            fusFs(results_gtype[,l,])
          })
          fu.fs.results <- do.call(rbind, fu.fs)

          #by population for each strataNames(results_gtype) and results_gtype[,,pops]
          fu.fs.pop <- lapply(strataNames(results_gtype), function(s){
            lapply(locNames(results_gtype), function(l){
              fusFs(results_gtype[,l,s])
            })
          })
          fu.fs.results.pop <- do.call(rbind, lapply(fu.fs.pop, function(x){
            do.call(rbind,x)
          }))

          fu.fs.all <- rbind(fu.fs.results, fu.fs.results.pop)

          # Tajimas D
          # by gene
          t.d <- lapply(locNames(results_gtype), function(l){
            tajimasD(results_gtype[,l,])
          })
          t.d.results <- do.call(c,t.d)

          # by gene per population
          t.d.pop <- lapply(strataNames(results_gtype), function(s){
            lapply(locNames(results_gtype), function(l){
              tajimasD(results_gtype[,l,s])
            })
          })
          t.d.pop.bind <- do.call(rbind, lapply(t.d.pop, function(x){
            do.call(rbind,x)
          }))
          t.d.all <- rbind(t.d.results, t.d.pop.bind)

          # Summary for loci and populations
          # Start Here - need summary for gene1 and gene2 not for pop1 and pop2??right
          smryLoci.gene <- lapply(locNames(results_gtype), function(l){
            summary(results_gtype[,l,])
          })

          smryLoci <- summary(results_gtype)
          smryLoci$strata.smry
          smryPop <- lapply(locNames(results_gtype), function(l){
            summary(results_gtype[,l,], by.strata = TRUE)$strata.smry
          })

          smryLP <- rbind(smryLoci$strata.smry, do.call(rbind, smryPop))

          # Nucleotide and percent within strata divergence
          # gene by population
          dA <- nucleotideDivergence(results_gtype)
          dA.names <- colnames(dA[[1]]$within)
          dA.pop <- do.call(rbind, lapply(1:length(dA), function(i){
            rbind(dA[[i]]$within)
          }))

          dA.all <- rbind(NA,NA, dA.pop)

          # More to do before can add
          # num.private.alleles
          hapFreqs <- lapply(locNames(results_gtype), function(l){
            hapFreqs <- alleleFreqs(results_gtype[,l,], by.strata = TRUE)
            by.loc <- sapply(hapFreqs, function(loc) {
              mat <- loc[, "freq", ]
              rowSums(apply(mat, 1, function(r) {
                result <- rep(FALSE, length(r))
                if(sum(r > 0) == 1) result[r > 0] <- TRUE
                result
                }))
              })
            colSums(by.loc)
          })
          hapFreqs.pop <- do.call(rbind,hapFreqs)

          #by gene
          hapFreqs.gene <- lapply(locNames(results_gtype), function(l){
            hapFreqs <- alleleFreqs(results_gtype[,l,], by.strata = FALSE)
            by.loc <- sapply(hapFreqs, function(loc) {
              mat <- loc[,"freq"] #no strata
              mat[mat>0]<-1
              sum(mat)
            })
            sum(by.loc)
          })
          hapFreqs.gene <- do.call(rbind,hapFreqs.gene)

          #Ne placeholder

          #how to deal with nD?
          results <- cbind(fu.fs.all,t.d.all,smryLP,dA.all)
          analysis_names <- c("fusFs",colnames(t.d.all),colnames(smryLP),dA.names,"num.private.haps")

          # Create the data array first time through
          if(is.null(params@analysis.results[["Locus"]][[curr_scn]])){
            params@analysis.results[["Locus"]][[curr_scn]] <- array(0, dim=c(num_loci*(num_pops+1),
                                                                             length(analysis_names),
                                                                             num_reps),
                                                                    dimnames = list(1:(num_loci*(num_pops+1)),
                                                                                    analysis_names,
                                                                                    1:num_reps))
          }

          params@analysis.results[["Locus"]][[curr_scn]][,,curr_rep] <-  locus.final
      }


    ###########################  Pairwise  ###########################

    if(params@analyses.requested["Pairwise"]){

      ############ START HERE ######################

      # genind

      if(inherits(params@rep.sample, "genind")){

        ##### no difference between pws and genotype data pws
        #Pairwise Chi2, D, F..., G...

        #testing nancycats - find and remove the NA locus
        sapply(locNames(results_gtype), function(l){
          results_gtype[,l,]
        })


        pairwiseTest(gene_gtype, nrep =5, keep.null=TRUE)
        # with nancycats example containing missing data: use poppr::missingno??
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
if(inherits(params@rep.sample,"multidna")){
  genes <- rep.sample$dna.seqs
  names(genes@dna) <- paste("gene", 1:length(genes@dna))
  id <- genes@labels
  df <- data.frame(id = id, strata = rep.sample[[1]], hap = id)
  rep.sample.gtypes <- df2gtypes(df, 1)

}

if(inherits(rep.sample,"DNAbin"))

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