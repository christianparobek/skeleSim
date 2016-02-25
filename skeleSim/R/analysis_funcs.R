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

  results_gtype <- params@rep.sample
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
  sA <- cbind(sA[, 1:2], mean = rowMeans(sA[, -(1:2), drop = FALSE], na.rm = TRUE), sA[, 3:ncol(sA)])
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
  cd <- if(ploidy(g) %in% 1:2) {
    dat <- genind2hierfstat(gtypes2genind(g))
    is.diploid <- ploidy(g) == 2
    chord.dist <- calcChordDist(dat, is.diploid)
    # chord.dist by locus
    chord.dist.locus <- lapply(locNames(g), function(l) {
      print(l)
      result <- calcChordDist(dat[, c("pop", l)], is.diploid)
      cbind(result[, 1:2], Locus = l, result[, 3])
    })
  } else NULL

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
