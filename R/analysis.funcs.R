#' @name analysis.funcs
#' @title Analysis functions for skeleSim parameter objects
#' @description Run Global, Locus, and Pairwise analyses on results from
#'   a single simulation replicate stored in params@rep.sample#'
#'
#' @param params a \linkS4class{skeleSim.params} object.
#' @param g a \linkS4class{gtypes} object.
#' @param num.perm.reps number of permutation replicates.
#' @param mat results matrix to be loaded into params object.
#' @param label analysis type label ("Global", "Locus", or "Pairwise").
#' @param dat data.frame in hierfstat format (see \code{\link[hierfstat]{genet.dist}}).
#'
#' @importFrom strataG statChi2 statFst statPhist statFstPrime statGst statGstPrime statGstDblPrime statFis
#'   overallTest pairwiseTest summarizeLoci hweTest mRatio
#'   fusFs nucleotideDivergence nucleotideDiversity privateAlleles
#'   sharedAlleles tajimasD theta strataSplit ploidy locNames
#'   strata<- gtypes2loci strataNames numAlleles
#' @importFrom reshape2 melt
#' @importFrom pegas hw.test
#'
#' @export
#'
analysisFunc <- function(params) {
  if(is.null(params@analysis.results)) {
    params@analysis.results <- lapply(1:length(params@scenarios), function(x) {
      list(Global = NULL, Locus = NULL, Pairwise = NULL)
    })
  }

  results.gtype <- params@rep.sample
# -->> REMOVE FOR RELEASE: SAVING gtypes OBJECT FOR TESTING <<--
  # label <- currentLabel(params)
  # file <- paste(label, ".results.gtype.rdata", sep = "")
  # save(results.gtype, file = file)
#-----
  num.perm.reps <- params@num.perm.reps

  opt <- options(warn = -1)

  if(params@analyses.requested["Global"]) {
    cat("  Global analysis...\n")
    mat <- globalAnalysis(results.gtype, num.perm.reps)
    params <- loadResultsMatrix(params, mat, "Global")
  }

  if(params@analyses.requested["Locus"]) {
    cat("  Locus analysis...\n")
    mat <- if(ploidy(results.gtype) > 1) {
      locusAnalysisGenotypes(results.gtype)
    } else {
      locusAnalysisHaplotypes(results.gtype)
    }
    params <- loadResultsMatrix(params, mat, "Locus")
  }

  if(params@analyses.requested["Pairwise"] & currentScenario(params)@num.pops > 1) {
    cat("  Pairwise analysis...\n")
    mat <- pairwiseAnalysis(results.gtype, num.perm.reps)
    params <- loadResultsMatrix(params, mat, "Pairwise")

  }

  options(opt)
  return(params)
}


#' @rdname analysis.funcs
#'
loadResultsMatrix <- function(params, mat, label) {
  curr_scn <- params@current.scenario
  num_reps <- params@num.sim.reps

  # alter mat to make sure all the rows are there.
  if (label=="Pairwise") {
    ###This first part is an overwrought way to specify all combinations of pops and loci
    np <- params@scenarios[[curr_scn]]@num.pops
    all.names.df <- expand.grid(
      pop1 = 1:np,
      pop2 = 1:np,
      locus = paste0("L", 1:params@scenarios[[curr_scn]]@num.loci)
    )
    all.names.df[, 1:2] <- t(apply(all.names.df[, 1:2], 1, sort))

    pops <- unique(all.names.df[, 1:2])
    pops <- pops[pops[, 1] != pops[, 2] ,]
    all.names.df <- merge(pops,all.names.df, all.x = T, all.y = F)
    all.names.df <- unique(all.names.df[with(all.names.df, order(pop1, pop2, locus)), ])
    ###this next line results in a complete list of pops and loci in the format that is used for rownames in the
    ###analysis.results pairwise matrices
    all.names <- switch(
      class(params@scenarios[[curr_scn]]@simulator.params),
      rmetasim.params = with(all.names.df, paste(pop1, pop2, locus, sep = "_")),
      fastsimcoal.params =  with(all.names.df, paste0("Sample ", pop1, "_Sample ", pop2, "_Locus_", gsub("L", "", locus)))
    )

    # eq zero means all predicted popxpopxloc are there
    all.present <- sum(sapply(all.names, function(x) {!x %in% rownames(mat)})) > 0
    if(all.present) {
      insrt <- all.names[which(!all.names %in% rownames(mat))]
      add <- matrix(NA, ncol = dim(mat)[2], nrow = length(insrt))
      rownames(add) <- insrt
      mat <- rbind(mat, add)
    }
    mat <- mat[order(rownames(mat)), ]
  }

  if(is.null(params@analysis.results[[curr_scn]][[label]])) {
    params@analysis.results[[curr_scn]][[label]] <- array(
      NA, dim = c(nrow(mat), ncol(mat), num_reps),
      dimnames = list(rownames(mat), colnames(mat), 1:num_reps)
    )
  }

  curr_rep <- params@current.replicate
  # print(colnames(params@analysis.results[[curr_scn]][[label]][, , curr_rep]))
  # print(colnames(mat))
  params@analysis.results[[curr_scn]][[label]][, , curr_rep] <- mat

  return(params)
}


#' @rdname analysis.funcs
#'
formatOverallStats <- function(g, num.perm.reps) {
  stats <- if(ploidy(g) == 1) {
    list(statChi2, statFst, statPhist)
  } else {
    list(statChi2, statFst, statFstPrime, statGst, statGstPrime,
         statGstDblPrime, statFis)
  }

  result <- overallTest(
    g, nrep = num.perm.reps, stats = stats, quietly = TRUE, max.cores = 1,
    model = "raw"
  )$result
  result.names <- paste(
    rep(rownames(result), each = 2), c("", ".pval"), sep = ""
  )
  result <- as.vector(t(result))
  names(result) <- result.names
  return(result)
}


#' @rdname analysis.funcs
#'
globalAnalysis <- function(g, num.perm.reps) {
  loc_names <- locNames(g)

  # run by locus analysis across all populations
  r.m <- lapply(loc_names, function(l) {
    formatOverallStats(g[, l, ], num.perm.reps)
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
  results.matrix <- rbind(
    formatOverallStats(g, num.perm.reps), results.matrix.l
  )
  analyses <- colnames(results.matrix)
  num_analyses <- length(analyses)
  rownames(results.matrix) <- c("Overall", loc_names)
  return(results.matrix)
}


#' @rdname analysis.funcs
#' @importFrom reshape2 melt
#' @importFrom pegas Fst
#'
locusAnalysisGenotypes <- function(g) {
  loc_names <- locNames(g)

  # by locus, all the other stats (num alleles etc) pulled from summarizeLoci
  smryLoci <- summarizeLoci(g)
  smryLoci <- data.frame(Pop = NA, Locus = rownames(smryLoci), smryLoci,
                         stringsAsFactors = FALSE)

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
  geno.cols <- which(!colnames(smry) %in% c("num.genotyped", "pct.genotyped"))
  smry <- smry[, geno.cols]
  rownames(smry) <- NULL


  save(file="g.rda",g)

  theta.hwe <- function(g) {
    cbind(theta = theta(g), hwe.p = hweTest(g, use.genepop = FALSE))
  }
  #cat("HWE tests started\n")
  th.locus <- data.frame(theta.hwe(g))
  #cat("theta.hwe run\n")
  th.locus <- cbind(
    Pop = NA, Locus = rownames(th.locus), th.locus, stringsAsFactors = FALSE
  )
  th.pop <- do.call(rbind, lapply(strataSplit(g), function(st.g) {
    df <- data.frame(theta.hwe(st.g))
    cbind(Pop = strataNames(st.g), Locus = rownames(df), df, stringsAsFactors = FALSE)
  }))
  th.all <- rbind(th.locus, th.pop)
  rownames(th.all) <- NULL

  #cat("HWE tests done\n")

  mratio.all <- NULL
  if(ploidy(g) == 2) {
    if(all(numAlleles(g) > 2)) { # diploid and not SNPs
      # mratio on gtypes object, function needs genetic data as a gtype
      mratio.locus <- mRatio(g, by.strata = FALSE, rpt.size = 1)
      mratio.all <- melt(t(mRatio(g, rpt.size = 1)))
      colnames(mratio.all) <- c("Pop", "Locus", "mRatio")
      mratio.all <- rbind(
        data.frame(
          Pop = NA, Locus = names(mratio.locus),
          mRatio = mratio.locus, stringsAsFactors = FALSE
        ),
        mratio.all
      )
      rownames(mratio.all) <- NULL
    } else {
      mratio.all <- rbind(
        data.frame(Pop = NA, Locus = locNames(g), mRatio = NA),
        expand.grid(Pop = strataNames(g), Locus = locNames(g), mRatio = NA)
      )
    }
  }

  # Number of private alleles by locus
  pa <- t(privateAlleles(g))
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
  FSTloci <- data.frame(
    Pop = NA, Locus = loc_names, FSTloci[loc_names, ], stringsAsFactors = FALSE
  )

  #for pops
  # pop.1/locus.1:num_loci - pop.num_pops/locus.1:num_loci..
  FSTpop <- lapply(levels(g.loci$population), function(x) {
    Fst(g.loci[g.loci$population == x,], pop = g.loci$population)
  })
  names(FSTpop) <- levels(g.loci$population)

  FSTpop <- do.call(rbind, mapply(function(mat, pop){
    data.frame(
      Pop = pop, Locus = loc_names, mat[loc_names, ], stringsAsFactors = FALSE
    )
  }, mat = FSTpop, pop = names(FSTpop), SIMPLIFY = FALSE))
  rownames(FSTpop) <- NULL

  FST.all <- rbind(FSTloci, FSTpop)
  rownames(FST.all) <- NULL

  #all the analyses get bound here
  # sorted by Loci across all populaitons,
  #   then Locus.1/Pop.1:Pop.num_pops ... Locus.num_loci/Pop.1:Pop.num_pops

  by.cols <- c("Pop", "Locus")
  smry <- merge(smry, th.all, by = by.cols, all = TRUE)
  if (!is.null(mratio.all))
      smry <- merge(smry, mratio.all, by = by.cols, all = TRUE)
  smry <- merge(smry, num.priv.allele, by = by.cols, all = TRUE)
  smry <- merge(smry, FST.all, by = by.cols, all = TRUE)
  smry$Pop <- as.character(smry$Pop)
  smry$Locus <- as.character(smry$Locus)
  smry <- smry[order(smry$Pop, smry$Locus), ]
  rownames(smry) <- sapply(1:nrow(smry), function(i) {
    if(is.na(smry$Pop[i])) {
      smry$Locus[i]
    } else {
      paste(smry$Locus[i], smry$Pop[i], sep = "_")
    }
  })
  smry$Pop <- smry$Locus <- NULL

  return(as.matrix(smry))
}


#' @rdname analysis.funcs
#'
hapSmryFunc <- function(g) {
  g <- g[, , , drop = TRUE]
  unstrat <- g
  strata(unstrat) <- "Default"
  t(sapply(locNames(unstrat), function(l) {
    loc.g <- unstrat[, l, , drop = TRUE]
    smry <- strataG::summary(loc.g)$strata.smry[1, ]
    smry <- smry[!names(smry) %in% "num.missing"]
    dvsty <- mean(nucleotideDiversity(g), na.rm = TRUE)
    Fs <- fusFs(loc.g)
    tD <- tajimasD(loc.g)[, "D"]
    c(smry, mean.nucleotide.diversity = dvsty, Fus.Fs = Fs, Tajimas.D = tD)
  }))
}


#' @rdname analysis.funcs
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


#' @rdname analysis.funcs
#' @importFrom hierfstat genet.dist
#' @importFrom utils combn
#'
calcChordDist <- function(dat) {
  chord.dist <- genet.dist(dat, diploid = TRUE, method = "Dch")
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


#' @rdname analysis.funcs
#' @importFrom hierfstat genind2hierfstat
#' @importFrom stats aggregate
#'
pairwiseAnalysis <- function(g, num.perm.reps) {
  # pairwise population structure test
  stats <- if(ploidy(g) == 1) {
    list(statChi2, statFst, statPhist)
  } else {
    list(statChi2, statFst, statFstPrime, statGst, statGstPrime,
         statGstDblPrime, statFis)
  }

  # print("in pairwise analysis")
  pws.all <- pairwiseTest(
    g, nrep = num.perm.reps, stats = stats, quietly = TRUE, max.cores = 1,
    model = "raw"
  )$result

  pws.all$pair.label <- pws.all$n.1 <- pws.all$n.2 <- NULL
  pws <- do.call(rbind, lapply(locNames(g), function(l) {
    result <- pairwiseTest(
      g[, l, ], nrep = num.perm.reps, stats = stats, quietly = TRUE,
      max.cores = 1, model = "raw"
    )$result
    result$pair.label <- result$n.1 <- result$n.2 <- NULL
    cbind(result[, 1:2], Locus = l, result[, 3:ncol(result), drop = FALSE])
  }))

  pws <- rbind(
    cbind(pws.all[, 1:2], Locus = NA, pws.all[, 3:ncol(pws.all), drop = FALSE]),
    pws
  )

  # shared alleles
  sA <- sharedAlleles(g)
  sA <- cbind(
    sA[, 1:2],
    mean = rowMeans(sA[, -(1:2), drop = FALSE], na.rm = TRUE),
    sA[, 3:ncol(sA), drop = FALSE]
  )
  sA <- melt(
    sA, id.vars = c("strata.1", "strata.2"),
    variable.name = "Locus", value.name = "shared.alleles"
  )
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
      cbind(dA.all[, 1:2], Locus = NA, dA.all[, 3:ncol(dA.all), drop = FALSE]),
      dA.locus
    )
  } else NULL

  # chord distance
  cd <- NULL
  # cd <- if(ploidy(g) == 2) {
  #   dat <- genind2hierfstat(gtypes2genind(g))
  #   chord.dist <- calcChordDist(dat)
  #   # chord.dist by locus
  #   chord.dist.locus <- do.call(rbind, lapply(colnames(dat)[-1], function(l) {
  #     if (length(unique(dat[, l])) > 1) {
  #       result <- calcChordDist(dat[, c("pop", l)])
  #     } else {
  #       result <- matrix(0, nrow = length(unique(dat[,"pop"])), ncol = 3)
  #       result <- as.data.frame(result)
  #       names(result) <- c("strata.1","strata.2","chord.distance")
  #     }
  #     cbind(result[, 1:2], Locus = l, result[, 3])
  #   }))
  #   colnames(chord.dist.locus)[4] <- "chord.dist"
  #   chord.dist.locus
  # } else NULL

  # combine results into single matrix
  by.cols <- c("strata.1", "strata.2", "Locus")
  smry <- merge(pws, sA, by = by.cols)
  if(!is.null(dA)) smry <- merge(smry, dA, by = by.cols)
  if(!is.null(cd)) smry <- merge(smry, cd, by = by.cols)
  smry <- smry[with(smry, order(Locus, strata.1, strata.2)), ]
  rownames(smry) <- sapply(1:nrow(smry), function(i) {
    st.pair <- paste(smry$strata.1[i], smry$strata.2[i], sep = "_")
    if(is.na(smry$Locus[i])) {
      st.pair
    } else {
      paste(st.pair, smry$Locus[i], sep = "_")
    }
  })
  smry$strata.1 <- smry$strata.2 <- smry$Locus <- NULL
  smry <- as.matrix(smry)

  return(smry)
}
