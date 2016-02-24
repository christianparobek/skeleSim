#' @name analysis_funcs
#' @title Analysis functions
#' @description Run Global, Locus, and Pairwise analyses on results from
#'   a single simulation replicate stored in params@rep.sample#'
#'
#' @param params a \linkS4class{skeleSim.params} object.
#' @param g a gtypes object.
#' @param num.perm.reps number of permutation replicates.
#' @param num.cores number of CPU cores to use.
#' @param mat results matrix to be loaded into params object.
#' @param label analysis type label ("Global", "Locus", or "Pairwise").
#' @param dat data.frame in hierfstat format (see \code{\link[hierfstat]{genet.dist}}).
#' @param is.diploid logical - is this a diploid object?
#'
#' @import strataG
#' @export
#'
analysis_funcs <- function(params){
  if(is.null(params@analysis.results)) {
    empty.list <- lapply(1:length(params@scenarios), function(x) NULL)
    params@analysis.results <- list(
      Global = empty.list, Locus = empty.list, Pairwise = empty.list
    )
  }

  results_gtype <- params@rep.results
  params@analyses.requested <- analyses.check(params@analyses.requested)
  num.perm.reps <- params@num.perm.reps
  num.cores <- params@num.cores

  opt <- options(warn = -1)

  if(params@analyses.requested["Global"]) {
    mat <- globalAnalysis(results_gtype, num.perm.reps, num.cores)
    params <- loadResultsMatrix(params, mat, "Global")
  }

  if(params@analyses.requested["Locus"]) {
    mat <- if(ploidy(results_gtype) > 1) {
      locusAnalysisGenotypes(results_gtype)
    } else {
      locusAnalysisHaplotypes(results_gtype)
    }
    params <- loadResultsMatrix(params, mat, "Locus")
  }

  if(params@analyses.requested["Pairwise"]) {
    mat <- pairwiseAnalysis(results_gtype, num.perm.reps, num.cores)
    params <- loadResultsMatrix(params, mat, "Pairwise")
  }

  ######################### Global ##########################

  if(params@analyses.requested["Global"]){

    # multiDNA - from apex package
    if(inherits(params@rep.sample,c("multidna","gtypes","list"))){

      #overall_stats() is from skeleSim.funcs.R
      #run by locus analysis across all populations
      # i.e. for each locus, run analysis across populations
      r.m <- lapply(locNames(results_gtype), function (l){
        gtypes_this_loc<-results_gtype[,l,]
        overall_stats(gtypes_this_loc)
      })
      #find complete list of column names
      #Just in case any didn't get computed above?
      analysis.names <- unique(do.call(c, lapply(r.m, function(x) names(x))))
      r.m <- lapply(r.m, function(x){
        missing <- setdiff(analysis.names, names(x))
        x[missing] <- NA
        names(x) <- analysis.names
        x
      })
      results.matrix.l <- do.call(rbind, r.m) # make that list into a matrix
      results.matrix <- rbind(overall_stats(results_gtype),results.matrix.l)
      analyses <- colnames(results.matrix)
      num_analyses <- length(analyses)
      #rownames(results.matrix) <- locNames(results_gtype)

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

    # if data came in genind form?
    if(inherits(params@rep.sample, "genind")){

      # Eric will improve overallTest {strataG} to deal with invariant loci
      # all statistics crash
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
  options(opt)
  return(params)
}

#' @rdname analysis_funcs
#'
loadResultsMatrix <- function(params, mat, label) {
  curr_scn <- params@current.scenario
  num_reps <- params@num.reps
  if(is.null(params@analysis.results[[label]][[curr_scn]])) {
    params@analysis.results[[label]][[curr_scn]] <- array(
      NA, dim = c(nrow(mat), ncol(mat), num_reps),
      dimnames = list(rownames(mat), colnames(mat), 1:num_reps)
    )
  }
  curr_rep <- params@current.replicate
  params@analysis.results[[label]][[curr_scn]][, , curr_rep] <- mat
  return(params)
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
      smryPop.1 <- cbind(locus_pop,do.call(rbind,summarizeLoci(results_gtype, by.strata = TRUE)))
      #sort by each population per locus
      smryPop <- smryPop.1[order(rownames(smryPop.1)),order(colnames(smryPop.1))]

      smry <- rbind(smryLoci[,-c(1:2)],smryPop[,-c(10,11)])

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
      FSTloci<-pegas::Fst(results_loci) ## had to add the pegas ref, because a circos error was in the way

      #for pops
      # pop.1/locus.1:num_loci - pop.num_pops/locus.1:num_loci..
      FSTpop<-lapply(levels(results_loci$population), function(x){
        fst <- pegas::Fst(results_loci[results_loci$population == x,], pop=results_loci$population)
      }) # had to add the pegas thing here too...
      FSTpop2sort <- mapply(function(mat,pn){
        x <- data.frame(mat)
        x$Locus <- row.names(x)
        x$Pop <- pn
        x
      }, mat = FSTpop, pn = 1:length(FSTpop), SIMPLIFY = FALSE)
      FSTpop.1 <- do.call(rbind,FSTpop2sort)
      FSTpop.2 <- FSTpop.1[order(FSTpop.1$Locus,FSTpop.1$Pop),]
      FSTpop.all <- rbind(FSTloci,FSTpop.2[,-c(4:5)])


      #all the analyses get bound here
      # sorted by Loci across all populaitons,
      #   then Locus.1/Pop.1:Pop.num_pops ... Locus.num_loci/Pop.1:Pop.num_pops
      locus.final <- cbind(HWE.pval = hwe[,-c(1:2)],mrat_results_all,smry,num.priv.allele,FSTpop.all, row.names=NULL)
      analysis_names <- colnames(locus.final)
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

#' @rdname analysis_funcs
#'
overall_stats <- function(g, num.perm.reps, num.cores) {
  ovl <- overallTest(g, nrep = num.perm.reps, num.cores = num.cores, quietly = TRUE)
  ovl.result <- ovl$result[complete.cases(ovl$result), ]
  global.wide <- as.vector(t(ovl.result))
  names(global.wide) <- paste(
    rep(rownames(ovl.result), each = 2), c("", ".pval"), sep = ""
  )
  return(global.wide)
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

#' @rdname analysis_funcs
#'
globalAnalysis <- function(g, num.perm.reps, num.cores) {
  loc_names <- locNames(g)

  # run by locus analysis across all populations
  r.m <- lapply(loc_names, function(l) {
    overall_stats(g[, l, ], num.perm.reps, num.cores)
  })
  # find complete list of column names
  analysis.names <- unique(unlist(lapply(r.m, names)))
  r.m <- lapply(r.m, function(x){
    missing <- setdiff(analysis.names, names(x))
    x[missing] <- NA
    names(x) <- analysis.names
    x
  })
  results.matrix.l <- do.call(rbind, r.m)
  results.matrix <- rbind(overall_stats(g, num.perm.reps, num.cores), results.matrix.l)
  analyses <- colnames(results.matrix)
  num_analyses <- length(analyses)
  rownames(results.matrix) <- c("Overall", loc_names)
  return(results.matrix)
}


#' @rdname analysis_funcs
#' @importFrom reshape2 melt
#' @importFrom pegas Fst
#'
locusAnalysisGenotypes <- function(g) {
  loc_names <- locNames(g)

  # by locus, all the other stats (num alleles etc) pulled from summarizeLoci
  smryLoci <- summarizeLoci(g)
  smryLoci <- data.frame(Pop = NA, Locus = rownames(smryLoci), smryLoci, stringsAsFactors = FALSE)

  # by population
  smryLociPop <- summarizeLoci(g, by.strata = TRUE)
  smryLociPop <- do.call(rbind, lapply(names(smryLociPop), function(pop) {
    result <- data.frame(
      Pop = pop, Locus = rownames(smryLociPop[[pop]]),
      smryLociPop[[pop]], stringsAsFactors = FALSE
    )
    rownames(result) <- NULL
    result
  }))

  smry <- rbind(smryLoci, smryLociPop)
  rownames(smry) <- NULL

  # mratio on gtypes object, function needs genetic data as a gtype
  mratio.all <- melt(t(mRatio(g)))
  colnames(mratio.all) <- c("Pop", "Locus", "mRatio")

  # Number of private alleles by locus
  pa <- privateAlleles(g)
  # this has the number of alleles that are private per locus
  perLocus <- colSums(pa)
  by.loc <- melt(pa)
  colnames(by.loc) <- c("Pop", "Locus", "num.priv.allele")
  perLocus <- data.frame(
    Pop = NA, Locus = names(perLocus), num.priv.allele = perLocus,
    stringsAsFactors = FALSE
  )
  # the rows will be have the private alleles for each population by locus
  num.priv.allele <- rbind(perLocus, by.loc)
  rownames(num.priv.allele) <- NULL

  # Convert from genind to loci (package pegas))
  #Fis estimation:   Fis = 1-(Ho/He)
  #for loci
  g.loci <- gtypes2loci(g)
  FSTloci <- Fst(g.loci)
  FSTloci <- data.frame(Pop = NA, Locus = loc_names, FSTloci[loc_names, ], stringsAsFactors = FALSE)

  #for pops
  # pop.1/locus.1:num_loci - pop.num_pops/locus.1:num_loci..
  FSTpop <- lapply(levels(g.loci$population), function(x) {
    Fst(g.loci[g.loci$population == x,], pop = g.loci$population)
  })
  names(FSTpop) <- levels(g.loci$population)

  FSTpop <- do.call(rbind, mapply(function(mat, pop){
    data.frame(Pop = pop, Locus = loc_names, mat[loc_names, ], stringsAsFactors = FALSE)
  }, mat = FSTpop, pop = names(FSTpop), SIMPLIFY = FALSE))
  rownames(FSTpop) <- NULL

  FST.all <- rbind(FSTloci, FSTpop)
  rownames(FST.all) <- NULL

  #all the analyses get bound here
  # sorted by Loci across all populaitons,
  #   then Locus.1/Pop.1:Pop.num_pops ... Locus.num_loci/Pop.1:Pop.num_pops

  locus.final <- merge(smry, mratio.all, by = c("Pop", "Locus"), all = TRUE)
  locus.final <- merge(locus.final, num.priv.allele, by = c("Pop", "Locus"), all = TRUE)
  locus.final <- merge(locus.final, FST.all, by = c("Pop", "Locus"), all = TRUE)
  locus.final$Pop <- as.character(locus.final$Pop)
  locus.final$Locus <- as.character(locus.final$Locus)
  locus.final <- locus.final[order(locus.final$Pop, locus.final$Locus), ]
  rownames(locus.final) <- sapply(1:nrow(locus.final), function(i) {
    if(is.na(locus.final$Pop[i])) {
      locus.final$Locus[i]
    } else {
      paste(locus.final$Locus[i], locus.final$Pop[i], sep = "_")
    }
  })
  locus.final$Pop <- locus.final$Locus <- NULL
  locus.final <- as.matrix(locus.final[order(rownames(locus.final)), ])

  return(locus.final)
}

    # If data comes in multidna form
    if(inherits(params@rep.sample, c("multidna","gtypes","list"))){


#' @rdname analysis_funcs
#'
hapSmryFunc <- function(g) {
  g <- g[, , , drop = TRUE]
  unstrat <- g
  strata(unstrat) <- "Default"
  smry <- t(sapply(locNames(g), function(l) {
    summary(unstrat[, l, , drop = TRUE])$strata.smry[1, ]
  }))
  dvsty <- sapply(nucleotideDiversity(g), mean, na.rm = TRUE)
  Fs <- fusFs(g)
  tD <- tajimasD(g)[, "D"]
  cbind(smry, mean.nucleotide.diversity = dvsty, Fus.Fs = Fs, Tajimas.D = tD)
}


#' @rdname analysis_funcs
#'
locusAnalysisHaplotypes <- function(g) {
  smry.all <- hapSmryFunc(g)
  smry.pop <- do.call(rbind, lapply(strataNames(g), function(st) {
    result <- hapSmryFunc(g[, ,st])
    rownames(result) <- paste(rownames(result), st, sep = "_")
    result
  }))

  # number of private alleles
  pa <- privateAlleles(g)
  pa.melt <- melt(pa)

  # Ne placeholder

  smry.all <- cbind(smry.all, num.private.alleles = rowSums(pa))
  smry.pop <- cbind(smry.pop, num.private.alleles = pa.melt[, 3])
  smry <- rbind(smry.all, smry.pop)

  return(smry)
}

    # If data comes in genind form
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

      locus.final <- cbind(pws.all[,-c(1:5)],sA = sA.all,chord_distance = chord.dist.all[,2])
      #locus.final <- locus.final.names[,sapply(locus.final.names,is.numeric)]
      analysis_names <- names(locus.final)

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

      params@analysis.results[["Pairwise"]][[curr_scn]][,,curr_rep] <-  as.matrix(as.matrix(locus.final))


#' @rdname analysis_funcs
#' @importFrom hierfstat genet.dist
#'
calcChordDist <- function(dat, is.diploid) {
  chord.dist <- genet.dist(dat, diploid = is.diploid, method = "Dch")
  chord.dist <- as.matrix(chord.dist)
  rownames(chord.dist) <- colnames(chord.dist) <- levels(dat$pop)
  pop.pairs <- data.frame(combn(levels(dat$pop), 2), stringsAsFactors = FALSE)
  result <- do.call(rbind, lapply(pop.pairs, function(pop.pair) {
    data.frame(
      strata.1 = pop.pair[1], strata.2 = pop.pair[2],
      chord.distance = chord.dist[pop.pair[1], pop.pair[2]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  return(result)
}


  }
  params # what does this do?


#' @rdname analysis_funcs
#'
pairwiseAnalysis <- function(g, num.perm.reps, num.cores) {
  # pairwise population structure test
  stats <- if(ploidy(g) == 1) {
    list(statChi2, statFst, statPhist)
  } else {
    list(statChi2, statFst, statFstPrime, statGst, statGstPrime,
         statGstDblPrime, statFis)
  }
  pws.all <- pairwiseTest(
    g, nrep = num.perm.reps, stats = stats,
    quietly = TRUE, num.cores = num.cores
  )$result
  pws.all$pair.label <- pws.all$n.1 <- pws.all$n.2 <- NULL
  pws <- do.call(rbind, lapply(locNames(g), function(l) {
    result <- pairwiseTest(
      g[, l, ], nrep = num.perm.reps, stats = stats,
      quietly = TRUE, num.cores = num.cores
    )$result
    result$pair.label <- result$n.1 <- result$n.2 <- NULL
    cbind(result[, 1:2], Locus = l, result[, 3:ncol(result)])
  }))
  pws <- rbind(cbind(pws.all[, 1:2], Locus = NA, pws.all[, 3:ncol(pws.all)]), pws)

  # shared alleles
  sA <- sharedAlleles(g)
  sA <- cbind(sA[, 1:2], mean = rowMeans(sA[, -(1:2)], na.rm = TRUE), sA[, 3:ncol(sA)])
  sA <- melt(sA, id.vars = c("strata.1", "strata.2"),
             variable.name = "Locus", value.name = "shared.alleles")
  sA$Locus[sA$Locus == "mean"] <- NA

  dA <- if(ploidy(g) == 1) {
    # nucleotide Divergence and mean.pct.between
    dA.locus <- lapply(nucleotideDivergence(g), function(x) x$between[, 1:4])
    dA.locus <- do.call(rbind, lapply(names(dA.locus), function(l) {
      x <- dA.locus[[l]]
      cbind(x[, 1:2], Locus = l, x[, 3:ncol(x)])
    }))
    colnames(dA.locus)[5] <- "mean.nucleotide.divergence"
    dA.all <- aggregate(
      dA.locus[, -(1:3)],
      list(strata.1 = dA.locus$strata.1, strata.2 = dA.locus$strata.2),
      mean, na.rm = TRUE
    )
    rbind(
      cbind(dA.all[, 1:2], Locus = NA, dA.all[, 3:ncol(dA.all)]), dA.locus
    )
  } else NULL

  #   # chord.dist
#   if(ploidy(g) %in% 1:2) {
#     dat <- genind2hierfstat(gtypes2genind(g))
#     is.diploid <- ploidy(g) == 2
#     chord.dist <- calcChordDist(dat, is.diploid)
#     # chord.dist by locus
#     chord.dist.locus <- lapply(locNames(g), function(l) {
#       print(l)
#       result <- calcChordDist(dat[, c("pop", l)], is.diploid)
#       cbind(result[, 1:2], Locus = l, result[, 3])
#     })
#   } else NULL

  smry <- merge(pws, sA, by = c("strata.1", "strata.2", "Locus"))
  if(!is.null(dA)) smry <- merge(smry, dA, by = c("strata.1", "strata.2", "Locus"))
  smry <- smry[with(smry, order(Locus, strata.1, strata.2)), ]
  rownames(smry) <- sapply(1:nrow(smry), function(i) {
    st.pair <- paste(smry$strata.1[i], smry$strata.2[i], sep = "_")
    if(is.na(smry$Locus[i])) st.pair else paste(st.pair, smry$Locus[i], sep = "_")
  })
  smry$strata.1 <- smry$strata.2 <- smry$Locus <- NULL
  smry <- as.matrix(smry)

  return(smry)
}
