set.commonparams<-function(main_list, gen_ind_list){
  hardcode=TRUE
  if(hardcode==TRUE){
    num_pops <- 4
    main_list$common_params$num_pops <- num_pops
    main_list$common_params$pop_sizes <- rep(1000, num_pops)
    main_list$common_params$overall_mig_rate <- 0.001
    spatially.explicit <- FALSE
    mig.model <- if(spatially.explicit) "distance" else "island"
    R.int <- if(spatially.explicit) {
      main_list$common_params$mean_mig_dist <- 2
      main_list$common_params$h.dim <- c(4, 5)
      if(prod(h.dim) != num_pops) stop("h.dim does is not compatible with num_pops")
      landscape.mig.matrix(h = num_pops,
                           h.dim = h.dim,
                           mig.model = mig.model,
                           distance.fun = dexp,
                           rate = 1 / main_list$common_params$mean_mig_dist,
                           distance.factor = 1
      )$R.int
    } else {
      landscape.mig.matrix(h = num_pops, mig.model = mig.model)$R.int
    }
    main_list$common_params$mig_rates <- R.int * main_list$common_params$overall_mig_rate
    main_list$common_params$marker_type<-factor("microsat",levels=c("microsat","snp","sequence"))
    main_list$common_params$num_markers<-3
    main_list$common_params$mut_rate<-0.0005
    main_list$common_params$sample_sizes<-rep(20,num_pops)

    return(main_list)
  }

}