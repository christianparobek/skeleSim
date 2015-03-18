set.commonparams<-function(main_list, gen_ind_list){
  hardcode=TRUE
  if(hardcode==TRUE){ 
    num_pops<-16
    main_list$common_params$num_pops<-num_pops
    main_list$common_params$pop_sizes<-rep(1000,num_pops)
    main_list$common_params$mig_model<-"distance"
    main_list$common_params$overall_mig_rate<-0.01
    main_list$common_params$mean_mig_dist<-2
    main_list$common_params$mig_rates<-landscape.mig.matrix(h=num_pops,h.dim=c(4,4),
            mig.model=main_list$common_params$mig_model,
            distance.fun=dexp,rate=1/main_list$common_params$mean_mig_dist,distance.factor=1
            )$R.int*main_list$common_params$overall_mig_rate
    main_list$common_params$marker_type<-factor("microsat",levels=c("microsat","snp","sequence"))
    main_list$common_params$num_markers<-10
    main_list$common_params$mut_rate<-rep(0.005,num_pops)
    main_list$common_params$sample_sizes<-rep(20,num_pops)
    
    return(main_list)
  }
  
}