#' @title Analysis functions
#' @description Run Global, Locus, and Pairwise analyses on results from
#'   a single simulation replicate stored in params@rep.sample#'
#'
#' @param params a \linkS4class{skeleSim.params} object.
#'
#' @export
#'
analysis_func <- function(params){
  if(is.null(params@analysis.results)) {
    params@analysis.results <- list(Global = NULL, Locus = NULL, Pairwise = NULL)
  }

  results_gtype <- results2gtypes(params)
  params@analyses.requested <- analyses.check(params@analyses.requested)

  if(params@analyses.requested["Global"]) {
    mat <- globalAnalysis(params, results_gtype)
    params < loadResultsMatrix(params, mat, "Global")
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
    mat <- pairwiseAnalysis(params, results_gtype)
    params <- loadResultsMatrix(params, mat, "Global")
  }

  return(params)
}


loadResultsMatrix <- function(params, mat, label) {
  curr_scn <- params@current_scenario
  num_reps <- params@num.reps
  if(is.null(params@analysis.results[[label]][[curr_scn]])) {
    empty.arr <- array(
      0, dim = c(nrow(mat), ncol(mat), num_reps),
      dimnames = list(rownames(mat), colnames(mat), 1:num_reps)
    )
    params@analysis.results[[label]][[curr_scn]] <- empty.arr
  }
  curr_rep <- params@current.replicate
  params@analysis.results[[label]][[curr_scn]][, , curr_rep] <- mat
  return(params)
}


overall_stats <- function(g) {
  opt <- options(warn = -1)
  ovl <- overallTest(g, nrep = 5, quietly = TRUE)
  ovl.result <- ovl$result[complete.cases(ovl$result), ]
  global.wide <- as.vector(t(ovl.result))
  names(global.wide) <- paste(
    rep(rownames(ovl.result), each = 2), c("", ".pval"), sep = ""
  )
  options(opt)
  global.wide
}


globalAnalysis <- function(params, g) {
  loc_names <- locNames(g)

  # run by locus analysis across all populations
  r.m <- lapply(loc_names, function(l) {
    overall_stats(g[, l, ])
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
  results.matrix <- rbind(overall_stats(g), results.matrix.l)
  analyses <- colnames(results.matrix)
  num_analyses <- length(analyses)
  rownames(results.matrix) <- c("Overall", loc_names)

  loadResultsMatrix(params, results.matrix, "Global")
}


#' @importFrom reshape2 melt
locusAnalysisGenotypes <- function(params, g) {
  opt <- options(warn = -1)

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
  perLocus <- colSums(by.loc)
  by.loc <- melt(by.loc)
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

  options(opt)
  return(locus.final)
}


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


locusAnalysisHaplotypes <- function(g) {
  opt <- options(warn = -1)

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

  options(opt)
  return(smry)
}


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
  result
}


pairwiseAnalysis <- function(g) {
  # nucleotide Divergence and mean.pct.between
  dA <- lapply(nucleotideDivergence(g), function(x) x$between[, 1:3])

  # pairwise test
  stats <- if(ploidy(g) == 1) {
    list(statChi2, statFst, statPhist)
  } else {
    list(statChi2, statFst, statFstPrime, statGst, statGstPrime,
         statGstDblPrime, statFis)
  }
  pws.all <- pairwiseTest(g, nrep = 5, stat.list = stats, quietly = TRUE)$result
  pws.all$pair.label <- pws.all$n.1 <- pws.all$n.2 <- NULL
  pws <- lapply(locNames(g), function(l) {
    result <- pairwiseTest(
      g[, l, ], nrep = 5, stat.list = stats, quietly = TRUE
    )$result
    result$pair.label <- result$n.1 <- result$n.2 <- NULL
    result
  })
  names(pws) <- locNames(g)

  # shared alleles
  sA <- sharedAlleles(g)
  sA.mean <- rowMeans(sA[, -(1:2)])
  sA.sum <- rowSums(sA[, -(1:2)])
  sA <- cbind(sA, mean = sA.mean, sum = sA.sum)
  sA <- melt(sA, id.vars = c("strata.1", "strata.2"),
       variable.name = "Locus", value.name = "shared.alleles")

  # chord.dist
  if(ploidy(g) %in% 1:2) {
    dat <- genind2hierfstat(gtypes2genind(g))
    is.diploid <- ploidy(g) == 2
    chord.dist <- calcChordDist(dat, is.diploid)
    # chord.dist by locus
    chord.dist.locus <- lapply(locNames(g), function(l) {
      print(l)
      result <- calcChordDist(dat[, c("pop", l)], is.diploid)
      cbind(result[, 1:2], Locus = l, result[, 3])
    })
  }
}
