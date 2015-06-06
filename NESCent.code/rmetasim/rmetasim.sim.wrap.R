#' @title rmetasim.sim.wrap
#' A wrapper for rmetasim simulations to be used in the skeleSim framework.
#' @parlist The main input parameter list for the simulation.  This must contain the
#' element common_params and the element spec_params_rmetasim
#'
rmetasim.sim.wrap <- function(parlist)
    {
        testsnp <- F
        testseq <- F
        
        ## this function assumes that loci are sequences of length 1.  A and C are lumped together and G and T
        ## are lumped together.  The allele calls are altered and the landscape$individual object is modified
        landscape.snp.convert <- function(land)
            {
                lvec <- landscape.locusvec(land)
                for (loc in 1:length(land$loci))
                    {
                        lcols <- which(lvec==loc)+6
                        states <-  landscape.locus.states(loc,land)
                        acind <- states$aindex[states$state %in% c("A","C")]
                        if (length(acind)>0)
                            {
                                alleles <- c(land$individuals[,lcols])
                                alleles[alleles %in% acind] <- acind[1]
                                land$individuals[,lcols] <- alleles
                            }
                        gtind <- states$aindex[states$state %in% c("G","T")]
                        if (length(gtind)>0)
                            {
                                alleles <- c(land$individuals[,lcols])
                                alleles[alleles %in% gtind] <- gtind[1]
                                land$individuals[,lcols] <- alleles
                            }
#                        if(!is.landscape(land)) stop("failed landscape test")
                    }
                land
            }

        ##do some error checking
        if (! is.list(parlist)) {stop("parlist must be a list object")}
        if (! ("common_params" %in% names(parlist))) {stop ("need a common_params list")}
        if (! ("spec_params_rmetasim" %in% names(parlist))) {stop ("need a spec_params_rmetasim list")}

        stgs <- dim(parlist$spec_params_rmetasim$surv_matrix)[1]
        habs <- parlist$common_params$num_pops
        l <- landscape.new.empty()
        l <- landscape.new.intparam(l,h=habs,s=stgs)
        l <- landscape.new.floatparam(l,s=parlist$spec_params_rmetasim$selfing_rate)
        l <- landscape.new.switchparam(l,mp=1)

        S <- parlist$spec_params_rmetasim$surv_matrix
        R <- parlist$spec_params_rmetasim$repr_matrix
        M <- parlist$spec_params_rmetasim$male_matrix
        l <- landscape.new.local.demo(l, S, R, M)

###set up landscape-level stuff
        carry <- parlist$common_params$pop_sizes
        extinct <- rep(0, length(carry))

        S <- matrix(0,habs*stgs,habs*stgs)
        R <- landscape.mig.matrix(h=habs,s=stgs,R.custom=parlist$common_params$mig_rates)$R
        M <- S
        l <- landscape.new.epoch(l, S=S,R=R,M=M,carry=carry,extinct=extinct)

        loctype <- as.character(parlist$common_params$locus_type)
if (testsnp)
    {
        loctype="snp"
        parlist$common_params$allele_size <- 1
        parlist$common_params$mut_rate <- rep(parlist$common_params$mut_rate,parlist$common_params$num_loci)
    }
if (testseq)
    {
        loctype="sequence"
        parlist$common_params$num_loci <- 1
        parlist$common_params$allele_size <- parlist$common_params$sequence_length
        parlist$common_params$mut_rate <- rep(parlist$common_params$mut_rate,parlist$common_params$num_loci)
    }
     

        
        ltype <- which(c("microsat","sequence","snp")==loctype)
        for (loc in 1:parlist$common_params$num_loci)
            {
                if (ltype==1) #ssr
                    {
                        l <- landscape.new.locus(l, type = 1, ploidy = 2, transmission=0,
                                                 mutationrate = parlist$common_params$mut_rate[loc],
                                                 numalleles = length(parlist$spec_params_rmetasim$init_allele_freqs[[loc]]),
                                                 frequencies = parlist$spec_params_rmetasim$init_allele_freqs[[loc]])
                    }
                else if (ltype==2) #seq
                    {
                        l <- landscape.new.locus(l,
                                                 type = 2, transmission = 1,
                                                 ploidy = 1, 
                                                 mutationrate = parlist$common_params$mut_rate[loc],
                                                 numalleles = length(parlist$spec_params_rmetasim$init_allele_freqs[[loc]]),
                                                 frequencies = parlist$spec_params_rmetasim$init_allele_freqs[[loc]],
                                                 allelesize = parlist$common_params$allele_size
                                                 )
                    }
                else if (ltype==3) #snp 
                    {
                        l <- landscape.new.locus(l,
                                                 type = 2, transmission=0,
                                                 ploidy = 2, 
                                                 mutationrate = parlist$common_params$mut_rate[loc],
                                                 numalleles = 4,
####                                                 frequencies = parlist$spec_params_rmetasim$init_allele_freqs[[loc]], 
                                                 allelesize = 1
                                                 )
                    }
            }
        
        l <- landscape.new.individuals(l,do.call(c,lapply(carry,function(x){x*rep(1/stgs,stgs)})))

        if (testsnp|testseq)
            {
                l <- landscape.simulate(l,20)
            } else {
                l <- landscape.simulate(l,parlist$spec_params_rmetasim$num_gens_to_run)
            }

        l <- landscape.sample(l,ns=parlist$common_params$sample_sizes)

        if (ltype==1)
            {
                parlist$rep.result <- landscape.make.genind(l)
            }
        else if (ltype==2)
            {
                states <- as.data.frame(landscape.locus.states(1,l))
                genos <- data.frame(pop=landscape.populations(l),aindex=l$individuals[,7])
                seq <- merge(genos,states,all.x=T)
                seq <- seq[order(seq$pop),]
                dna.seq <- strsplit(as.character(tolower(seq$state)),"")
                parlist$rep.result <- list(strata=data.frame(seq$pop),
                                           dna.seq=as.DNAbin(do.call(rbind,strsplit(tolower(as.character(seq$state)),""))))
            }
        else if (ltype==3)
            {
                parlist$rep.result <- landscape.make.genind(landscape.snp.convert(l))
            }
        parlist
    }

