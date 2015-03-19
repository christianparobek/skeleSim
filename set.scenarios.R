##########################
#  Set scenario function #
###########################

#this function will make a grid of scenarios to try 
#this grid will be all the multipliers to try, e.g.
#to try twice the population sizes

set.scenarios<-function(main_list, samp_size_mult, num_loc_mult, pop_size_mult, mig_rate_mult, mut_rate_mult, ...){
  #hardcoding in values for these multipliers
  #later they will be input by the user
  samp_size_mult<-c(2,5)
  num_loc_mult<-1
  pop_size_mult<-1
  mig_rate_mult<-c(1,2)
  mut_rate_mult<-1  
  main_list$scenarios_list<-
              expand.grid("samps_mult"=samp_size_mult,"loc_mult"=num_loc_mult,"pop_size_mult"=pop_size_mult,
              "mig_mult"=mig_rate_mult, "mut_mult"=mut_rate_mult)
  
  return(main_list)
  
}