#' @title is.megalist
#' @megalist a megalist object
#' Performs some internal megalist checks.  Returns a boolean and 
#'
is.megalist <- function(megalist)
    {
        correct <- TRUE #flag for correct format and contents

        megalist_elements <- c("proj_title", "proj_date" , "goal_chosen","sim_chosen", "user_has_data", "common_params", "analysis_func", "sim_func", "scenarios_list")
        common_names <- c("num_pops", "sample_sizes", "num_loci", "pop_sizes", "overall_mig_rate", "mig_rates", "locus_type", "mut_rate", "sequence_length","num_reps", "current_scenario", "current_replicate")
        
        cp <- megalist$common
        
        if (correct)
            {
                if (!is.list(megalist))
                   {
                       message("megalist should be a list")
                       correct <- FALSE
                   }
            }

        if (correct)
            {
                missing.names <- megalist_elements[!megalist_elements%in% names(megalist)]
                if (length(missing.names)>0)
                    {
                        message (paste("These elements are missing from the megalist",paste(missing.names,collapse=",")))
                        correct <- FALSE
                    }
            }
        
        if (correct)
            {
                missing.names <- common_names[!common_names %in% names(cp)]
                if (length(missing.names)>0)
                    {
                        message (paste("These names are missing from the common_param component of the megalist",paste(missing.names,collapse=",")))
                        correct <- FALSE
                    }
            }

        if (correct)
            {
                if (!is.function(megalist$analysis_func))
                    {
                        message("megalist$analysis_func must be a function")
                    }
            }
        
        if (correct)
            {
                if (!is.function(megalist$sim_func))
                    {
                        message("megalist$sim_func must be a function")
                    }
            }
        if (correct)
            {
                if (!(is.data.frame(megalist$scenarios_list)|is.matrix(megalist$scenarios)))
                    {
                        message("megalist$scenarios_list should be a matrix or data frame")
                        correct <-  FALSE
                    }
                if ((correct) & (!is.null(rownames(megalist$scenarios_list))))
                    {
                        message("megalist$scenarios_list should be a matrix or data frame")
                        correct <-  FALSE
                    }
            }

### now that the simple tests have been conducted, let's perform some checks on the parameter values
        #check that the correct specific params section is included based on other params
        if (correct) 
            {
                if (megalist$sim_chosen=="c")
                    {
                        if (!"spec_params_fastsimcoal" %in% names(megalist))
                            {
                                message("you chose a coalescent simulator but there are not any params for that type of simulator")
                                correct <- FALSE
                            }
                        else
                            {
                                message("you chose a forward time simulator but there are not any params for that type of simulator")
                                correct <- FALSE
                            }
                    }
            }
        correct
    }
