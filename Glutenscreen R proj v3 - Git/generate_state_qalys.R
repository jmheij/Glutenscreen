#project: Glutenscreen markov
#action: function to generate qalys for generate_net_benefit function
#title: generate_state_qalys.R




generate_state_qalys <- function(input_parameters) {
  n_samples <- dim(input_parameters)[1]
  
  starting_age <- 2
  n_cycles <- 82 - starting_age

  # define the QALYS associated with the states per cycle
  # There is one for each PSA sample and each state
  
  state_qalys <- array(dim=c(n_samples, n_cycles, n_states), dimnames = list(NULL, NULL, state_names))
  
  starting_age_columnandrow <- read.csv("brismod/data/starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  starting_age_row <- which(starting_age_columnandrow$Starting.age  > starting_age)[1] - 1
  
  eq5d_norms <- read.csv("eq5d_norms_nl.csv")
  eq5d_norms$age <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  
  
  n_agecategories <- ceiling((n_cycles/10) - 1)
  for(i_age_category in c(0:1)) { ##note that 0:1 does not match exactly with ages 2-18 but is very close
    age_category_indices <- (c(1:10) + i_age_category * 10)
    age_category_indices <- age_category_indices[age_category_indices <= n_cycles]

    state_qalys[, age_category_indices, "diagnosed no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_GFD_children
    state_qalys[, age_category_indices, "Undiagnosed no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_undiagnosedCD_children - (input_parameters$probability_late_diagnosis_children * input_parameters$probability_biopsy_children * input_parameters$disutility_biopsy_children)
  }
  for(i_age_category in c(2:n_agecategories)) { 
    age_category_indices <- (c(1:10) + i_age_category * 10)
    age_category_indices <- age_category_indices[age_category_indices <= n_cycles]
    
    state_qalys[, age_category_indices, "diagnosed no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_GFD_adults
    state_qalys[, age_category_indices, "Undiagnosed no complications"] <- eq5d_norms[starting_age_row + i_age_category, 2] * input_parameters$utility_undiagnosedCD_adults - (input_parameters$probability_late_diagnosis_adults * input_parameters$probability_biopsy_adults * input_parameters$disutility_biopsy_adults)
  }
  
  state_qalys[, , "diagnosed osteoporosis"] <- state_qalys[, , "diagnosed no complications"] - input_parameters$disutility_osteoporosis 
  state_qalys[, , "diagnosed NHL"] <- state_qalys[, , "diagnosed no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "diagnosed GIC"] <- state_qalys[, , "diagnosed no complications"] - input_parameters$disutility_GIC
  state_qalys[, , "Undiagnosed osteoporosis"] <- state_qalys[, , "Undiagnosed no complications"] - input_parameters$disutility_osteoporosis
  state_qalys[, , "Undiagnosed NHL"] <- state_qalys[, , "Undiagnosed no complications"] - input_parameters$disutility_NHL
  state_qalys[, , "Undiagnosed GIC"] <- state_qalys[, , "Undiagnosed no complications"] - input_parameters$disutility_GIC
  state_qalys[, , "Death"] <- 0
  return(state_qalys[, , ])
}