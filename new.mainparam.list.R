new.mainparam.list<-function(title){
	list_of_params<-list(
		proj_title=as.character(title),
		goal_chosen=factor(x=NULL, levels=c("power","performance","nulldistribution")),
		sim_chosen=NULL,
		common_params=NULL,
		spec_params_fastsimcoal=NULL,
		spec_params_rmetasim=NULL	
		)
	return(list_of_params)

}