set.specparams.rmetasim<-function(main_list, gen_ind_list){
  #hardcode means these values are set- later we will allow them to be set by the user via interface
  hardcode=TRUE
  if(hardcode==TRUE){ 
    #extinction rate is vector of length number populations
    main_list$spec_params_rmetasim$ext_rates<-rep(0,main_list$common_params$num_pops)
    #hardcoded number of stages as 2, e.g. annual plant
    num_stages<-2
    #hardcorded maternal output number of offspring
    offsp_per_fem<-10
    #what follows are repr, surv and male repr matrices for rmetasim
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
    main_list$spec_params_rmetasim$num_gens_to_run<-100
    main_list$spec_params_rmetasim$max_pop_size<-20000
    main_list$spec_params_rmetasim$selfing_rate<-0.0
    
    
  }
}