#' @title rmetasim.sim.wrap
#' A wrapper for rmetasim simulations to be used in the skeleSim framework.
#' @parlist The main input parameter list for the simulation.  This must contain the
#' element common_params and the element spec_params_rmetasim
#'
rmetasim.sim.wrap <- function(parlist)
    {
        #function to convert rmetasim to genind object (takes care of sequences)xs
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
        loctype="snp"
        ltype <- which(c("microsat","sequence","snp")==loctype)
        for (loc in 1:parlist$common_params$num_loci)
            {
                if (ltype==1)
                    {
                        l <- landscape.new.locus(l, type = 1, ploidy = 2, 
                                                 mutationrate = parlist$common_params$mut_rate[loc],
                                                 numalleles = length(parlist$spec_params_rmetasim$init_allele_freqs[[loc]]),
                                                 frequencies = parlist$spec_params_rmetasim$init_allele_freqs[[loc]])
                    } else {
                        l <- landscape.new.locus(l,
                                                 type = 2,
                                                 ploidy = 2, 
                                                 mutationrate = parlist$common_params$mut_rate[loc],
                                                 numalleles = length(parlist$spec_params_rmetasim$init_allele_freqs[[loc]]),
                                                 frequencies = parlist$spec_params_rmetasim$init_allele_freqs[[loc]],
                                                 allelesize = parlist$common_params$allele_size
                                                 )
                    }
            }
        l <- landscape.new.individuals(l,do.call(c,lapply(carry,function(x){x*rep(1/stgs,stgs)})))
        l <- landscape.simulate(l,parlist$spec_params_rmetasim$num_gens_to_run)
        if (ltype==1)
            {
                landscape.make.genind(l)
            }
        else if (ltype==2)
            {
                NULL  #make a DNAbin object
            }
        else if (ltype==3)
            {
                NULL #make a SNP genind
            }
        
    }

