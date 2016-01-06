function(params){

  stopifnot(require(strataG))
  stopifnot(require(pegas))
  stopifnot(require(hierfstat))

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


  loc_names <- locNames(results_gtype)
  strata_names <- strataNames(results_gtype)

  # If analysis results is empty, the first analysis done creates the list to hold the data
  if(is.null(params@analysis.results)){
    params@analysis.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)
  }

        ######################### Global ##########################

  # TO DO: Need an across all loci row if a genind object in addition to each loci Will need a different
  # analysis.results format

  if(params@analyses.requested["Global"]){

    # multiDNA
    if(inherits(params@rep.sample,c("multidna","gtypes"))){

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
      params@analysis.results[["Global"]][[curr_scn]] <- array(0,dim=c(num_loci+1,
                                                                       num_analyses,
                                                                       num_reps),
                                                               dimnames = list(c(1:num_loci+1),
                                                                               analyses,
                                                                               1:num_reps))
    }

    params@analysis.results[["Global"]][[curr_scn]][,,curr_rep] <- results.matrix
    }

    if(inherits(params@rep.sample, "genind")){

      ovl <- lapply(locNames(results_gtype), function(l){
        overallTest(results_gtype[,l,], nrep = 5, stat.list = statList("chi2"), quietly = TRUE)$result
      })

      ovl.all <- lapply(ovl, function(x){
        as.data.frame(as.table(x))
        })

      ovl.all.values <- lapply(ovl.all, function(x){
        x[,-c(1:2)]
        })
      ovl.all.names <- paste(ovl.all[[1]]$Var1, ovl.all[[1]]$Var2,sep="_")
      ovl.all.out <- do.call(rbind, ovl.all.values)

      if(is.null(params@analysis.results[["Global"]][[curr_scn]])){
        params@analysis.results[["Global"]][[curr_scn]] <- array(0,dim=c(num_loci+1,
                                                                         num_analyses,
                                                                         num_reps),
                                                                 dimnames = list(c(1:num_loci+1),
                                                                                 analyses,
                                                                                 1:num_reps))
      }

      params@analysis.results[["Global"]][[curr_scn]][,,curr_rep] <- results.matrix
    }

    }


        ######################### LOCUS ############################

  if(params@analyses.requested["Locus"]){

# genind
        if(inherits(params@rep.sample,"genind")){

          #Hardy Weinberg, per locus over populations, per locus per population
          strataSplit(results_gtype)
          locus <- hweTest(results_gtype)
          locus.pop <- lapply(strataSplit(results_gtype), function(s){
            hw <- hweTest(s)
            missingloc <- setdiff(loc_names, names(hw))
            hw[missingloc] <- NA
            hw
            })
          hwe.pval <- do.call(rbind,lapply(locus.pop, function(x) x[match(names(locus.pop[[1]]), names(x))]))


          #mratio on gtypes object, function needs genetic data as a gtype
          mrat_results_all<-calc.mratio(results_gtype)
          matrix(mrat_results_all,num_pops,num_loci)

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


          # Convert from genind to loci (package pegas))
          #Fis estimation:   Fis = 1-(Ho/He)
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
            mean(nucleotideDiversity(results_gtype[,l,]@sequences))
          })
          nD <- do.call(rbind, r.m.gene)


          # by gene per popualation "strata"
          r.m <- lapply(strataNames(results_gtype), function(s){
            lapply(locNames(results_gtype), function(l){
            mean(nucleotideDiversity(results_gtype[,l,s]@sequences))
            })
          })
          r.m.bind <- do.call(c, lapply(r.m, function(x){
            do.call(c,x)
          }))
          nD.all <- c(nD, r.m.bind)

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
              if(is.null(fusFs(results_gtype[,l,s]))){
                NA
              } else {
              fusFs(results_gtype[,l,s])
            }
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
          t.d.results <- do.call(rbind,t.d)

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
          # need strata.smry data, not sequence summary...
          smryLoci.gene <- lapply(locNames(results_gtype), function(l){
            summary(results_gtype[,l,])#$seq.smry
          })

          smryLoci <- do.call(rbind,smryLoci.gene)

          smryPop <- lapply(locNames(results_gtype), function(l){
            summary(results_gtype[,l,], by.strata = TRUE)$strata.smry
          })
          smryPop.all <- do.call(rbind, smryPop)
          summary.analyses <- dimnames(smryPop.all)[[2]]
          #Loci over all populations, locus 1 per population, locus 2 per population...
          smryLP <- rbind(smryLoci, do.call(rbind, smryPop))

          # Nucleotide and percent within strata divergence
          # gene by population
          dA <- nucleotideDivergence(results_gtype)
          dA.names <- colnames(dA[[1]]$within)
          dA.pop <- do.call(rbind, lapply(1:length(dA), function(i){
            rbind(dA[[i]]$within)
          }))

          geneNAs <- matrix(NA, num_loci, 6)  # no data over strata for each gene
          dA.all <- rbind(geneNAs, dA.pop)


          # num.private.alleles
          hapFreqs <- lapply(strataNames(results_gtype), function(s){
            lapply(locNames(results_gtype), function(l){
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
          num.pri.haps <- rbind(hapFreqs.gene, hapFreqs.pop)

          #Ne placeholder

          #how to deal with nD?
          # make sure nucleotide Divergence is right and add or keep names below
          # to do start here
          results <- cbind(fu.fs.all,t.d.all,smryLP,dA.all)
          analysis_names <- c("nucloetide.diversity", "Fu.F",colnames(t.d.all),summary.analyses,
                              "nucleotide.divergence","mean.pct.within","num.private.haps", "Ne")

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
  }


        ######################### Pairwise ##########################

  if(params@analyses.requested["Pairwise"]){

      # genind

      if(inherits(params@rep.sample, "genind")){

        ##### no difference between pws and genotype data pws
        #Pairwise Chi2, D, F..., G...
        split.foo <- sapply(locNames(results_gtype), function(l){
          if(summary(results_gtype[,l,])$locus.smry[3] == 0){
            1
          } else {
            0
          }
        })

        summary(results_gtype[,4,])$allele.freqs
        summary(results_gtype[,4,])$locus.smry

        results_gtype_narm <- results_gtype[,locNames(results_gtype) != "fca45",]
        ############ START HERE ######################
        # fca45 has NaN for nancycats
        # using microbov excample

        pws <- pairwiseTest(results_gtype, nrep =5, keep.null=TRUE, quietly = TRUE)
        pws[1]$result #PHist + p.vals always NA
        sts <- c("pair.label","Chi2","Chi2.p.val","Fst","Fst.p.val","PHIst","PHIst.p.val")
        pws.sts <- pws[1]$result[sts]
        pws.sts$pair.label <- gsub("\\s*\\([^\\)]+\\)","",as.character(pws.sts$pair.label))
        pws.sts$strata.1 <- gsub(" .*$", "", as.character(pws.sts$pair.label))
        pws.sts$strata.2 <- gsub(".*v.","", as.character(pws.sts$pair.label))

        #pairwise order differs from pws.sts
        sA <- sharedAlleles(results_gtype)

        #chord.dist
        results_hierfstat <- genind2hierfstat(params@rep.sample)
        chord.dist <- genet.dist(results_hierfstat, diploid = TRUE, method = "Dch")
        nrow(as.data.frame(as.table(chord.dist)))

        psw.out <- cbind(pws.sts[with(pws.sts, order(strata.1,strata.2)),],
                         sA[with(sA,order(strata.1,strata.2)),])
        psw.out[,grep("strata", names(psw.out), invert=TRUE)]

        analysis_names <- names(psw.out[-1])

      #Data.frame of summary data into simulation replicate
        # Create the data array first time through
        if(is.null(params@analysis.results[["Pairwise"]][[curr_scn]])){
          params@analysis.results[["Pairwise"]][[curr_scn]] <- array(0, dim=c(length(data.frame(combn(1:num_pops,2))),
                                                                           length(analysis_names),
                                                                           num_reps),
                                                                  dimnames = list(1:length(data.frame(combn(1:num_pops,2))),
                                                                                  analysis_names,
                                                                                  1:num_reps))
        }

        params@analysis.results[["Locus"]][[curr_scn]][,,curr_rep] <-  locus.final

      }


    if(inherits(params@rep.sample, c("multidna","gtypes"))){

      #nucleotide Divergence and mean.pct.between
      dA <- nucleotideDivergence(results_gtype)
      dA.between <- lapply(dA, function(x){
        x$between
      })
      dA.all <- do.call(rbind,dA.between)

      #Chi2, Fst, PHist
      psw <- pairwiseTest(results_gtype,nrep = 5,stat.list = list(statGst),
                          quietly = TRUE)$result[-c(1:6,8:9,12:21)]

      # shared haps
      sA <- sharedAlleles(results_gtype)


      }
  params
    }
  }


