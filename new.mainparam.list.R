new.mainparam.list<-function(title="test default"){
	list_of_params<-list(
		proj_title=as.character(title),
    proj_date=date(),
		goal_chosen=factor(x=NULL, levels=c("power","performance","nulldistribution")),
		sim_chosen=NULL,
    user_has_data=FALSE,
		common_params=NULL,
		spec_params_fastsimcoal=NULL,
		spec_params_rmetasim=NULL,
    analyses_to_run=NULL,
    results_from_analyses=NULL
		)
	return(list_of_params)

}