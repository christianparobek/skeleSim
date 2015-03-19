set.commonparams<-function(main_list, gen_ind_obj=NULL){

    main_list$user_has_data <- !is.null(gen_ind_obj)

    hardcode <- TRUE
    if(hardcode){
        num_pops<-3
                                        #If the user has data (in genind format)...
                                        #these are parameters we can pull from their genind object metadata
                                        #right now it calls the metadata.getter function multiple times
                                        #could be more efficient to call once... EFFICIENCY
        if (main_list$user_has_data) {
            main_list$common_params$num_pops<-genind.metadata.getter(gen_ind_obj)$NumberOfPops
            main_list$common_params$sample_sizes<-genind.metadata.getter(gen_ind_obj)$SampsPerPop
            main_list$common_params$num_loci<-genind.metadata.getter(gen_ind_obj)$NumberOfLoci
        } else {
            main_list$common_params$num_pops<-num_pops
            main_list$common_params$sample_sizes<-rep(20,num_pops)
            main_list$common_params$num_loci<-10
        }
                                        #size of each population
        main_list$common_params$pop_sizes<-rep(1000,num_pops)

                                        # migration model for rmetasim
        spatially.explicit <- FALSE
        main_list$common_params$overall_mig_rate <- 0.01
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

        main_list$common_params$locus_type<-factor("microsat",levels=c("microsat","snp","sequence"))
        main_list$common_params$mut_rate<-0.0005
        main_list$common_params$sequence_length <- 400
                                        # number of simulation replicates to run
        main_list$common_params$num_reps<-100
        main_list$common_params$current_scenario<-1
        main_list$common_params$current_replicate<-1

        return(main_list)
}
