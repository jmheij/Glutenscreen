#project: Glutenscreen markov
#action: function to generate costs for generate_net_benefit function
#title: generate_state_costs.R

generate_state_costs_b <- function(input_parameters) {

  state_costs <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  starting_age_columnandrow <- read.csv("brismod/data/starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  starting_age_column <- which(starting_age_columnandrow$Starting.age  > starting_age)[1] - 1
  
  
  #copy probability_IDA for inclusion in state costs
  probability_IDA <- data.frame(input_parameters$probability_IDA_0, input_parameters$probability_IDA_10, input_parameters$probability_IDA_20, input_parameters$probability_IDA_30, input_parameters$probability_IDA_40, 
                                input_parameters$probability_IDA_50, input_parameters$probability_IDA_60, input_parameters$probability_IDA_70, input_parameters$probability_IDA_80, input_parameters$probability_IDA_90)
  
 
  n_agecategories <- ceiling((n_cycles/10) - 1)
  
  for(i_age_category in c(0:1)) {
    age_category_indices <- (c(1:10) + i_age_category * 10)
    age_category_indices <- age_category_indices[age_category_indices <= n_cycles]
    
    state_costs[, age_category_indices, "diagnosed no complications"] <- input_parameters$cost_CDGFD_children + 
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) + 
      input_parameters$cost_gfp 
    
    state_costs[, age_category_indices, "diagnosed osteoporosis"]  <- input_parameters$cost_osteoporosis + 
      input_parameters$cost_CDGFD_children + 
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) +
      input_parameters$cost_gfp 
    
    state_costs[, age_category_indices, "diagnosed GIC"]  <- input_parameters$cost_GIC + 
      input_parameters$cost_CDGFD_children +
      input_parameters$cost_gfp +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)

    state_costs[, age_category_indices, "diagnosed NHL"] <-  input_parameters$cost_NHL +
      input_parameters$cost_CDGFD_children + 
      input_parameters$cost_gfp +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)
    state_costs[, age_category_indices, "Undiagnosed no complications"] <- input_parameters$cost_undiagnosedCD_children + 
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) +
      (input_parameters$probability_late_diagnosis_children * input_parameters$cost_diagnosis_ns) 
    state_costs[, age_category_indices, "Undiagnosed osteoporosis"] <- input_parameters$cost_undiagnosedCD_children + 
      input_parameters$cost_osteoporosis + 
      (input_parameters$probability_late_diagnosis_children * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
    state_costs[, age_category_indices, "Undiagnosed GIC"] <- input_parameters$cost_GIC + 
      input_parameters$cost_undiagnosedCD_children + 
      (input_parameters$probability_late_diagnosis_children * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
   
    state_costs[, age_category_indices, "Undiagnosed NHL"] <- input_parameters$cost_NHL +
      input_parameters$cost_undiagnosedCD_children + 
      (input_parameters$probability_late_diagnosis_children * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)
  }
  
  for(i_age_category in c(2:n_agecategories)) {
    age_category_indices <- (c(1:10) + i_age_category * 10)
    age_category_indices <- age_category_indices[age_category_indices <= n_cycles]
    
    state_costs[, age_category_indices, "diagnosed no complications"] <- input_parameters$cost_CDGFD_adults + 
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) +
      input_parameters$cost_gfp 
    
    state_costs[, age_category_indices, "diagnosed osteoporosis"]  <- input_parameters$cost_osteoporosis + 
      input_parameters$cost_CDGFD_adults + 
      input_parameters$cost_gfp +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
    
    state_costs[, age_category_indices, "diagnosed GIC"]  <- input_parameters$cost_GIC + 
      input_parameters$cost_CDGFD_adults +
      input_parameters$cost_gfp +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)

    state_costs[, age_category_indices, "diagnosed NHL"] <-  input_parameters$cost_NHL +
      input_parameters$cost_CDGFD_adults + 
      input_parameters$cost_gfp +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)
    state_costs[, age_category_indices, "Undiagnosed no complications"] <- input_parameters$cost_undiagnosedCD_adults +
      (input_parameters$probability_late_diagnosis_adults * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
    
    state_costs[, age_category_indices, "Undiagnosed osteoporosis"] <- input_parameters$cost_undiagnosedCD_adults + 
      input_parameters$cost_osteoporosis + 
      (input_parameters$probability_late_diagnosis_adults * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 
    
    state_costs[, age_category_indices, "Undiagnosed GIC"] <- input_parameters$cost_GIC + 
      input_parameters$cost_undiagnosedCD_adults + 
      (input_parameters$probability_late_diagnosis_adults * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA) 

    state_costs[, age_category_indices, "Undiagnosed NHL"] <- input_parameters$cost_NHL +
      input_parameters$cost_undiagnosedCD_adults + 
      (input_parameters$probability_late_diagnosis_adults * input_parameters$cost_diagnosis_ns) +
      (probability_IDA[, starting_age_column + i_age_category] * input_parameters$cost_IDA)
  }
  
  
  state_costs[, , "Death"] <- 0
  return(state_costs[, , ])
}

