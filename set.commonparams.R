set.commonparams<-function(main_list, gen_ind_obj){
  hardcode=TRUE
  if(hardcode==TRUE){ 
    num_pops<-16`
    #If the user has data (in genind format)...
    #these are parameters we can pull from their genind object metadata
    #right now it calls the metadata.getter function multiple times
    #could be more efficient to call once... EFFICIENCY
    if (main_list$user_has_data==TRUE) {
      main_list$common_params$num_pops<-genind.metadata.getter(gen_ind_obj)$NumberOfPops
      main_list$common_params$sample_sizes<-genind.metadata.getter(gen_ind_obj)$SampsPerPop
      main_list$common_params$num_loci<-genind.metadata.getter(gen_ind_obj)$NumberOfLoci
    }
    if (main_list$user_has_data==FALSE) {
      main_list$common_params$num_pops<-num_pops
      main_list$common_params$sample_sizes<-rep(20,num_pops)
      main_list$common_params$num_loci<-10
    }
    #total number of replicate simulations to run
    main_list$common_params$pop_sizes<-rep(1000,num_pops)
    #migration model for rmetasim
    main_list$common_params$mig_model<-"distance" #choices are...
    #base or overall migration rate
    main_list$common_params$overall_mig_rate<-0.01
    main_list$common_params$mean_mig_dist<-2
    #ERROR CHECKING GOES HERE for h.dim/ num_pops... island, linear, distance, etc
    main_list$common_params$mig_rates<-landscape.mig.matrix(h=num_pops,h.dim=c(4,4),
            mig.model=main_list$common_params$mig_model,
            distance.fun=dexp,rate=1/main_list$common_params$mean_mig_dist,distance.factor=1
            )$R.int*main_list$common_params$overall_mig_rate
    main_list$common_params$locus_type<-factor("microsat",levels=c("microsat","snp","sequence"))
    main_list$common_params$mut_rate<-rep(0.005,num_pops)
    main_list$common_params$num_reps<-100
    main_list$common_params$current_scenario<-1
    
    return(main_list)
  }
  
}