#' @title Analysis functions
#' @description Analyses for genind or multidna data
#'
#' @param params Take results of a simulation from params@rep.sample and run Global, Locus, and Pairwise analyses
#'
#' @import strataG
#' @import pegas
#' @import hierfstat
#'
#' @export
analysis_funcs <- function(params){

  # saving global variables
  curr_scn<-params@current.scenario
  curr_rep<-params@current.replicate
  num_loci<-params@scenarios[[curr_scn]]@num.loci
  num_reps<-params@num.reps
  num_pops<-params@scenarios[[curr_scn]]@num.pops
  params@analyses.requested <- analyses.check(params@analyses.requested)

  #params@rep.sample is either a genind or a list of DNAbin objects
  if(inherits(params@rep.sample, "genind")){
    results_gtype<-genind2gtypes(params@rep.sample)
  } else if (inherits(params@rep.sample, "gtypes")){
    results_gtype <- params@rep.sample
  } else {
    # Convert the list of DNAbin objects to gtypes
    genes <- params@rep.sample
    results_gtype <- sequence2gtypes(genes$dna.seqs, strata = genes$strata)
    results_gtype <- labelHaplotypes(results_gtype)$gtype
  }

  loc_names <- locNames(results_gtype)
  strata_names <- strataNames(results_gtype)

  # If analysis results is empty, the first analysis done creates the list to hold the data
  if(is.null(params@analysis.results)){
    params@analysis.results <- sapply(c("Global","Locus","Pairwise"), function(x) NULL)
  }

  ######################### Global ##########################

  if(params@analyses.requested["Global"]){

    # multiDNA
    if(inherits(params@rep.sample,c("multidna","gtypes","list"))){

      #overall_stats() is from skeleSim.funcs.R
      #run by locus analysis across all populations
      r.m <- lapply(locNames(results_gtype), function (l){
        gtypes_this_loc<-results_gtype[,l,]
        overall_stats(gtypes_this_loc)
      })
      #find complete list of column names
      analysis.names <- unique(do.call(c, lapply(r.m, function(x) names(x))))
      r.m <- lapply(r.m, function(x){
        missing <- setdiff(analysis.names, names(x))
        x[missing] <- NA
        names(x) <- analysis.names
        x
      })
      results.matrix.l <- do.call(rbind, r.m)
      results.matrix <- rbind(overall_stats(results_gtype),results.matrix.l)
      analyses <- colnames(results.matrix)
      num_analyses <- length(analyses)
      rownames(results.matrix) <- c("OverLoci",loc_names)

      # We are printing by gene, not overall gene analysis. This differs from the genind code below.
      # The first row will hold summary statistics over all loci regardless of population structure.
      # The remaining rows will hold summary statistics per locus
      if(is.null(params@analysis.results[["Global"]][[curr_scn]])){
        params@analysis.results[["Global"]][[curr_scn]] <- array(0,dim=c(num_loci+1,
                                                                         num_analyses,
                                                                         num_reps),
                                                                 dimnames = list(c(1:(num_loci+1)),
                                                                                 analyses,
                                                                                 1:num_reps))
      }

      params@analysis.results[["Global"]][[curr_scn]][,,curr_rep] <- results.matrix
    }

    if(inherits(params@rep.sample, "genind")){

      ovl.global <- overallTest(results_gtype, nrep=5, stat.list=statList("chi2"), quietly=TRUE)$result
      ovl.all.global <- as.data.frame(as.table(ovl.global))[,3]

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
      ovl.all.results <- rbind(ovl.all.global,ovl.all.out)
      row.names(ovl.all.results) <- c("overall",loc_names)

      if(is.null(params@analysis.results[["Global"]][[curr_scn]])){
        params@analysis.results[["Global"]][[curr_scn]] <- array(0,dim=c(num_loci+1,
                                                                         length(ovl.all.names),
                                                                         num_reps),
                                                                 dimnames = list(c(1:(num_loci+1)),
                                                                                 ovl.all.names,
                                                                                 1:num_reps))
      }

      params@analysis.results[["Global"]][[curr_scn]][,,curr_rep] <- ovl.all.results
    }

  }


  ######################### LOCUS ############################

  if(params@analyses.requested["Locus"]){

    # genind
    if(inherits(params@rep.sample,"genind")){

      #Hardy Weinberg, per locus over populations, per locus per population
      locus <- hweTest(results_gtype)
      locus.pop <- lapply(strataSplit(results_gtype), function(s){
        hw <- hweTest(s)
        missingloc <- setdiff(loc_names, names(hw))
        hw[missingloc] <- NA
        hw
      })
      hwe.pval <- do.call(rbind,lapply(locus.pop, function(x) x[match(names(locus.pop[[1]]), names(x))]))
      hwe.pval.pop <- as.data.frame(as.table(hwe.pval))
      locus.t <- cbind(NA, names(locus), data.frame(locus))
      names(locus.t) <- c("Pop","Locus","HWE.pval")
      names(hwe.pval.pop) <- c("Pop","Locus","HWE.pval")
      hwe <- rbind(locus.t,hwe.pval.pop)

      #mratio on gtypes object, function needs genetic data as a gtype
      mrat_results_all<-calc.mratio(results_gtype)

      #by locus, all the other stats (num alleles etc) pulled from summarizeLoci
      smryLoci <- cbind(Locus = 1:num_loci,Pop = NA,summarizeLoci(results_gtype))
      #by population
      locus_pop <- as.matrix(expand.grid(1:num_loci,1:num_pops))  #cycles through the loci for each pop
      smryPop.1 <- data.frame(locus_pop,do.call(rbind,summarizeLoci(results_gtype, by.strata = TRUE)))
      #sort by each population per locus
      smryPop <- smryPop.1[with(smryPop.1, order(smryPop.1$Var1,smryPop.1$Var2)),]

      smry <- rbind(smryLoci[,-c(1:2)],data.matrix(smryPop[,-c(1,2)]))


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
      num.priv.allele <- c(perLocus, as.data.frame(as.table(by.loc))[,3])

      # Convert from genind to loci (package pegas))
      #Fis estimation:   Fis = 1-(Ho/He)
      results_loci<-genind2loci(params@rep.sample)
      #for loci
      FSTloci<-Fst(results_loci)

      #for pops
      # pop.1/locus.1:num_loci - pop.num_pops/locus.1:num_loci..
      FSTpop<-lapply(levels(results_loci$population), function(x){
        fst <- Fst(results_loci[results_loci$population == x,], pop=results_loci$population)
      })
      FSTpop2sort <- mapply(function(mat,pn){
        x <- data.frame(mat)
        x$Locus <- row.names(x)
        x$Pop <- pn
        x
      }, mat = FSTpop, pn = 1:length(FSTpop), SIMPLIFY = FALSE)
      FSTpop.1 <- do.call(rbind,FSTpop2sort)
      FSTpop.1 <- data.frame(locus_pop,FSTpop.1)
      FSTpop.2 <- FSTpop.1[order(FSTpop.1$Var1,FSTpop.1$Var2),]
      FSTpop.all <- rbind(FSTloci,FSTpop.2[,grep("F",colnames(FSTpop.2),value=TRUE)])


      #all the analyses get bound here
      # sorted by Loci across all populaitons,
      #   then Locus.1/Pop.1:Pop.num_pops ... Locus.num_loci/Pop.1:Pop.num_pops
      locus.final <- data.frame(HWE.pval = hwe[,-c(1:2)],mrat_results_all,smry,num.priv.allele,FSTpop.all, row.names=NULL)
      analysis_names <- colnames(locus.final)
      row.names(locus.final) <- c(loc_names,apply(expand.grid(1:num_pops,1:num_loci),1,
                                                  function(x) paste(x[2],x[1],sep="_")))
      locus.final <- as.matrix(locus.final)


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
    if(inherits(params@rep.sample,c("multidna","gtypes","list"))){

      #Nucleotide diversity
      # by gene, across populations
      r.m.gene <- lapply(locNames(results_gtype), function(l){
        mean(nucleotideDiversity(results_gtype[,l,]@sequences),na.rm=TRUE)
      })
      nD <- do.call(rbind, r.m.gene)

      # by gene per popualation "strata" - pop1:gene1, pop1:gene2, pop2....
      r.m <- lapply(strataNames(results_gtype), function(s){
        lapply(locNames(results_gtype), function(l){
          mean(nucleotideDiversity(results_gtype[,l,s]@sequences),na.rm=TRUE)
        })
      })
      r.m.bind <- do.call(c, lapply(r.m, function(x){
        do.call(c,x)
      }))
      nD.all <- c(nD, r.m.bind)

      # Fu's Fs
      fu.fs.results <- fusFs(results_gtype)

      #by population for each strataNames(results_gtype) and results_gtype[,,pops]
      fu.fs.pop <- lapply(strataNames(results_gtype), function(s){
        lapply(locNames(results_gtype), function(l){
          fusFs(results_gtype[,l,s])
        })
      })
      fu.fs.results.pop <- do.call(rbind, lapply(fu.fs.pop, function(x){
        do.call(rbind,x)
      }))

      fu.fs.all <- c(fu.fs.results, fu.fs.results.pop)

      # Tajimas D
      # by gene
      t.d.results <-  tajimasD(params@rep.sample$dna.seqs)

      # Num samples, num missing, num alleles, percent unique alleles, heterozygosity
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
      unstrat <- results_gtype
      strata(unstrat) <- "Default"
      smryLoci.gene <- lapply(locNames(unstrat), function(l){
        summary(unstrat[,l,])$strata.smry
      })
      smryLoci <- do.call(rbind,smryLoci.gene)

      smryPop <- lapply(locNames(results_gtype), function(l){
        summary(results_gtype[,l,], by.strata = TRUE)$strata.smry
      })
      smryPop.all <- do.call(rbind, smryPop)
      summary.analyses <- dimnames(smryPop.all)[[2]]

      #Loci over all populations, locus 1 per population, locus 2 per population...
      smryLP <- rbind(smryLoci,smryPop.all)

      # Nucleotide and percent within strata divergence, mean percent within
      # gene by population
      #mean.pct.within
      #          dA <- nucleotideDivergence(results_gtype)
      #          dA.names <- colnames(dA[[1]]$within)
      #          dA.pop <- do.call(rbind, lapply(1:length(dA), function(i){
      #            rbind(dA[[i]]$within)
      #          }))

      # do we want dA in locus, doesn't make sense.
      #          colnames(dA[[2]]$between)

      # nucleotide divergence is pairwise between and within strata. Cannot use unstratified.
      #          dA.genes <- nucleotideDivergence(unstrat)[[1]]$within
      #          geneNAs <- matrix(NA, num_loci, 6)  # no data over strata for each gene
      #          dA.all <- rbind(geneNAs, dA.pop)

      #          dA.all <- rbind(geneNAs, dA.pop)
      #          dA.analyses <- dimnames(dA.all)[[2]]


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
      hapFreqs.pop <- do.call(c,do.call(c,hapFreqs))

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
      num.pri.haps <- c(do.call(c,hapFreqs.gene), hapFreqs.pop)

      #Ne placeholder

      #how to deal with nD?
      # make sure nucleotide Divergence is right and add or keep names below
      # to do start here
      locus.final <- cbind(nD.all,fu.fs.all,t.d.all,smryLP,#dA.all[,1],
                           num.pri.haps)
      analysis_names <- c("nucloetide.diversity", "Fu.F",colnames(t.d.all),summary.analyses,
                          #"nucleotide.divergence",
                          "num.private.haps")
      row.names(locus.final) <- c(loc_names,
                                  apply(expand.grid(loc_names,strata_names),1,
                                        function(x) paste(x[2],x[1],sep="_")))

      # Create the data array first time through
      # gene.1...gene.num_loci, pop.1/gene.1:gene.num_loci...pop.num.pops/gene.1:gene.num_loci
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

    if(inherits(params@rep.sample, c("multidna","gtypes","list"))){

      #nucleotide Divergence and mean.pct.between
      dA <- nucleotideDivergence(results_gtype)
      dA.between <- lapply(dA, function(x){
        x$between
      })
      dA.all <- do.call(rbind,dA.between)[,grep("pct",names(dA.between[[1]]),invert = TRUE)]

      #Chi2, Fst, PHist
      psw <- lapply(locNames(results_gtype), function(l){
        pairwiseTest(results_gtype[,l,],nrep = 5,
                     stat.list = list(statGst),
                     quietly = TRUE)$result[-c(1:6,8:9,12:21)]
      })
      psw.all <- do.call(rbind,psw)

      # shared haps
      sA <- sharedAlleles(results_gtype)
      sA <- reshape(sA, direction = "long", varying=list(loc_names), v.names="sA",
                    timevar = "gene",idvar=c("strata.1","strata.2"),
                    new.row.names = 1:(num_loci*choose(num_pops,2)))[,"sA"]

      pairwise.final <- cbind(dA.all,psw.all,sA)
      analysis_names <- names(pairwise.final)[-c(1:2)]
      row.names(pairwise.final) <- apply(expand.grid(c(apply(combn(strata_names,2),2,
                                           function(x) paste(x[1],x[2],sep="_"))),loc_names),1,
                                         function(x) paste(x[2],x[1],sep=""))

      #Data.frame of summary data into simulation replicate
      # Create the data array first time through
      if(is.null(params@analysis.results[["Pairwise"]][[curr_scn]])){
        params@analysis.results[["Pairwise"]][[curr_scn]] <- array(0, dim=c(num_loci*choose(num_pops,2),
                                                                            length(analysis_names),
                                                                            num_reps),
                                                                   dimnames = list(apply(pairwise.final[,1:2],1,paste,collapse="."),
                                                                                   analysis_names,
                                                                                   1:num_reps))
      }

      params@analysis.results[["Pairwise"]][[curr_scn]][,,curr_rep] <-  data.matrix(pairwise.final[,-c(1:2)])

    }

    if(inherits(params@rep.sample, "genind")){

      #Pairwise Chi2, D, F..., G...
      pws <- pairwiseTest(results_gtype, nrep =5, keep.null=TRUE, quietly = TRUE)[1]$result
      #pairwise by locus
      pws.loc <- lapply(loc_names, function(l){
        pairwiseTest(results_gtype[,l,], nrep = 5, keep.null=TRUE, quietly=TRUE)[1]$result
      })
      pws.loc.all <- do.call(rbind,pws.loc)
      pws.all <- rbind(pws,pws.loc.all)

      sA <- sharedAlleles(results_gtype)
      sA.long <- reshape(sA, idvar = c("strata.1","strata.2"),
                         varying = names(sA[,-c(1:2)]),
                         timevar = "Locus",
                         v.names = "sA",
                         direction = "long")
      #shared alleles sumed over loci
      sA.sum <- rowSums(sA[,-c(1:2)])
      sA.all <- c(sA.sum,sA.long$sA)

      #chord.dist
      results_hierfstat <- genind2hierfstat(params@rep.sample)
      chord.dist <- genet.dist(results_hierfstat, diploid = TRUE, method = "Dch")
      ch.dist <- as.data.frame(as.table(chord.dist))
      # chord.dist by locus
      chord.dist.locus <- lapply(loc_names, function(l){
        as.data.frame(as.table(genet.dist(results_hierfstat[,c("pop",l)], diploid = TRUE, method = "Dch")))
      })
      chord.dist.l <- do.call(rbind,chord.dist.locus)
      chord.dist.all <- rbind(ch.dist,chord.dist.l)

      pws.final <- cbind(pws.all[,-c(1:5)],sA = sA.all,chord_distance = chord.dist.all[,2])
      #locus.final <- locus.final.names[,sapply(locus.final.names,is.numeric)]
      analysis_names <- names(pws.final)
      #### almost!!! ####
      row.names(pws.final) <- mapply(function(loc,pop){paste(loc,pop,sep=":")},
                                     loc_names,apply(combn(1:4,2),2,function(x) paste(x[1],x[2],sep="_")))

      #Data.frame of summary data into simulation replicate
      # Create the data array first time through
      if(is.null(params@analysis.results[["Pairwise"]][[curr_scn]])){
        params@analysis.results[["Pairwise"]][[curr_scn]] <- array(0, dim=c((num_loci*choose(num_pops,2))+choose(num_pops,2),
                                                                            length(analysis_names),
                                                                            num_reps),
                                                                   dimnames = list(1:(num_loci*choose(num_pops,2)+choose(num_pops,2)),
                                                                                   analysis_names,
                                                                                   1:num_reps))
      }

      params@analysis.results[["Pairwise"]][[curr_scn]][,,curr_rep] <-  as.matrix(pws.final)

    }


  }
  params
}


