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
  if(params@analyses.requested["Global"]) params <- globalAnalysis(params, results_gtype)
  if(params@analyses.requested["Locus"]) {
    params <- if(ploidy(results_gtype) > 1) {
      locusAnalysisGenotypes(params, results_gtype)
    } else {
      locusAnalysisHaplotypes(params, results_gtype)
    }
  }
  if(params@analyses.requested["Pairwise"]) params <- pairwiseAnalysis(params, results_gtype)

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
  params@analysis.results[[label]][[curr_scn]][, , params@current.replicate] <- mat
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
  af <- alleleFreqs(g, by.strata = TRUE)
  by.loc <- sapply(af, function(loc) {
    mat <- loc[, "freq", ]
    rowSums(apply(mat, 1, function(r) {
      result <- rep(FALSE, length(r))
      if(sum(r > 0) == 1) result[r > 0] <- TRUE
      result
    }))
  })
  rownames(by.loc) <- strataNames(g)
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

  loadResultsMatrix(params, locus.final, "Locus")
}

locusAnalysisHaplotypes(params, g) {
  # Nucleotide diversity by gene, across populations
  mean.nD <- do.call(rbind, lapply(locNames(g), function(l) {
    mean(nucleotideDiversity(g[, l, ]), na.rm = TRUE)
  }))

  # by gene per popualation "strata" - pop1:gene1, pop1:gene2, pop2....
  mean.nD.all <- do.call(c, sapply(stratSplit(g), function(st){
    do.call(c, lapply(locNames(st), function(l) {
      mean(nucleotideDiversity(st[, l, ]), na.rm = TRUE)
    }))
  }))

  # Fu's Fs
  Fs.results <- fusFs(g)

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

  loadResultsMatrix(params, locus.final, "Locus")
}


pairwiseAnalysis <- function(params, g) {
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
      #### check order of 3+ populations ####
      row.names(pws.final) <- c(apply(combn(1:num_pops,2),2,function(x){
        paste(x[1],x[2],sep="_")
      }),apply(expand.grid(loc_names,
                           apply(combn(1:num_pops,2),2,function(x) paste(x[1],x[2],sep="_"))),
               1,paste,collapse="."))

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
