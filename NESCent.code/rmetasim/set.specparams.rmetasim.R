set.specparams.rmetasim<-function(main_list, gen_ind_obj=NULL){
  
  #hardcode means these values are set- later we will allow them to be set by the user via interface
  hardcode=TRUE
  
  if(hardcode==TRUE){ 
    #extinction rate is vector of length number populations
    main_list$spec_params_rmetasim$ext_rates<-rep(0,main_list$common_params$num_pops)
    #hardcoded number of stages as 2, e.g. annual plant
    num_stages<-2
    #hardcorded maternal output/ number of offspring per female
    offsp_per_fem<-10
    #what follows are female reproduction, survival and male reproductive matrices for rmetasim
    #for either 'annual' or 'perennial' organisms
    if (num_stages==2){
      main_list$spec_params_rmetasim$repr_matrix<-matrix(c(0.0,offsp_per_fem,
                                                           0.0,0.0),nrow=2,byrow=T)
      main_list$spec_params_rmetasim$surv_matrix<-matrix(c(0.00,0.0,
                                                           0.5,0.0),nrow=2,byrow=T)
      main_list$spec_params_rmetasim$male_matrix<-matrix(c(0.0,0.0,
                                                            0.0,1.0),nrow=2,byrow=T)
    }
    if (num_stages==3){
      main_list$spec_params_rmetasim$repr_matrix<-matrix(c(0.0,1.0,offsp_per_fem,
                                                           0.0,0.0,0.0,
                                                           0.0,0.0,0.0),nrow=3,byrow=T)
      main_list$spec_params_rmetasim$surv_repr_matrix<-matrix(c(0.02,0.0,0.0,
                                                           0.2,0.6,0.0,
                                                           0.0,0.2,0.8),nrow=3,byrow=T)
      main_list$spec_params_rmetasim$male_matrix<-matrix(c(0.0,0.0,0.0,
                                                            0.0,1.0,1.0,
                                                            0.0,1.0,1.0),nrow=3,byrow=T)
    }
    #total number of generations to run
    main_list$spec_params_rmetasim$num_gens_to_run<-100
    #maximum rmetasim landscape size
    main_list$spec_params_rmetasim$max_pop_size<-20000
    #rmetasim selfing rate
    main_list$spec_params_rmetasim$selfing_rate<-0.0
    #initial allele frequencies for starting population
    if (is.null(gen_ind_obj))
        {
                main_list$spec_params_rmetasim$init_allele_freqs<-lapply(1:main_list$common_params$num_loci,
                                                                         function(loc)
                                                                         {
                                                                             nall <- rpois(1,lambda=3)+1
                                                                             frq <- runif(nall)
                                                                             frq <- frq/sum(frq)
                                                                             frq
                                                                         })
        } else {
                main_list$spec_params_rmetasim$init_allele_freqs<-genind.metadata.getter(gen_ind_obj)$FreqByLocus
        }

    
    return(main_list)
  }
}
