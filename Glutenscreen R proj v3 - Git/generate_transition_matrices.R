
#project: Glutenscreen markov
#action: function to generate transition matrices
#title: generate_transition_matrices.R


generate_transition_matrices <- function(input_parameters) {
  
  n_samples <- dim(input_parameters)[1]
  
  transition_matrices <- array(dim = c(n_samples, n_cycles, n_states, n_states),
                               dimnames = list(NULL, NULL, state_names, state_names))
  
  # There is one transition matrix for each cycle and each PSA sample
  # Store them in an array with (before filling in below) NA entries
  # Define transition matrices
  
  starting_age_columnandrow <- read.csv("brismod/data/starting_age_column.csv")
  starting_age_columnandrow$starting_age_row <- c(1:11)
  #Initial cohort at diagnosis - depends on age at diagnosis
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_column <- which(starting_age_columnandrow$Starting.age  > starting_age)[1] - 1
  
  
  
  # Osteoporosis probabilities 
  # Generate probabilities based on these rates and ratios
  osteoporosis_probability_GFD <- osteoporosis_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    osteoporosis_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("osteoporosis_lograte_", 10*(i_age_category - 1))] + input_parameters$log_or_osteoporosis_GFD))
    osteoporosis_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("osteoporosis_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_or_osteoporosis_noGFD))
  }
  osteoporosis_probability_GFD_all <- as.data.frame(osteoporosis_probability_GFD)
  colnames(osteoporosis_probability_GFD_all) <- c("osteoporosis_probability_GFD_0","osteoporosis_probability_GFD_10","osteoporosis_probability_GFD_20","osteoporosis_probability_GFD_30",
                                                  "osteoporosis_probability_GFD_40","osteoporosis_probability_GFD_50","osteoporosis_probability_GFD_60",
                                                  "osteoporosis_probability_GFD_70","osteoporosis_probability_GFD_80","osteoporosis_probability_GFD_90")
  osteoporosis_probability_noGFD_all <- as.data.frame(osteoporosis_probability_noGFD)
  colnames(osteoporosis_probability_noGFD_all) <- c("osteoporosis_probability_noGFD_0","osteoporosis_probability_noGFD_10","osteoporosis_probability_noGFD_20","osteoporosis_probability_noGFD_30",
                                                    "osteoporosis_probability_noGFD_40","osteoporosis_probability_noGFD_50","osteoporosis_probability_noGFD_60",
                                                    "osteoporosis_probability_noGFD_70","osteoporosis_probability_noGFD_80","osteoporosis_probability_noGFD_90")
  
  
  # Calculate NHL probabilities
  
  NHL_probability_GFD <- NHL_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    NHL_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("NHL_lograte_", 10*(i_age_category - 1))] + input_parameters$log_rr_NHL_GFD))
    NHL_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("NHL_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_rr_NHL_noGFD))
  }
  NHL_probability_GFD_all <- as.data.frame(NHL_probability_GFD)
  colnames(NHL_probability_GFD_all) <- c("NHL_probability_GFD_0","NHL_probability_GFD_10","NHL_probability_GFD_20","NHL_probability_GFD_30",
                                                  "NHL_probability_GFD_40","NHL_probability_GFD_50","NHL_probability_GFD_60",
                                                  "NHL_probability_GFD_70","NHL_probability_GFD_80","NHL_probability_GFD_90")
  NHL_probability_noGFD_all <- as.data.frame(NHL_probability_noGFD)
  colnames(NHL_probability_noGFD_all) <- c("NHL_probability_noGFD_0","NHL_probability_noGFD_10","NHL_probability_noGFD_20","NHL_probability_noGFD_30",
                                                    "NHL_probability_noGFD_40","NHL_probability_noGFD_50","NHL_probability_noGFD_60",
                                                    "NHL_probability_noGFD_70","NHL_probability_noGFD_80","NHL_probability_noGFD_90")
  
  #calulate GIC probabilities
  GIC_probability_GFD <- GIC_probability_noGFD <- matrix(NA, nrow = n_samples, ncol = 10)
  
  for(i_age_category in 1:10) {
    GIC_probability_GFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("GIC_lograte_", 10*(i_age_category - 1))] + input_parameters$log_rr_GIC_GFD))
    GIC_probability_noGFD[, i_age_category] <- 1 - exp(-exp(input_parameters[, paste0("GIC_lograte_", 10 * (i_age_category - 1))] + input_parameters$log_rr_GIC_noGFD))
  }
  GIC_probability_GFD_all <- as.data.frame(GIC_probability_GFD)
  colnames(GIC_probability_GFD_all) <- c("GIC_probability_GFD_0","GIC_probability_GFD_10","GIC_probability_GFD_20","GIC_probability_GFD_30",
                                         "GIC_probability_GFD_40","GIC_probability_GFD_50","GIC_probability_GFD_60",
                                         "GIC_probability_GFD_70","GIC_probability_GFD_80","GIC_probability_GFD_90")
  GIC_probability_noGFD_all <- as.data.frame(GIC_probability_noGFD)
  colnames(GIC_probability_noGFD_all) <- c("GIC_probability_noGFD_0","GIC_probability_noGFD_10","GIC_probability_noGFD_20","GIC_probability_noGFD_30",
                                           "GIC_probability_noGFD_40","GIC_probability_noGFD_50","GIC_probability_noGFD_60",
                                           "GIC_probability_noGFD_70","GIC_probability_noGFD_80","GIC_probability_noGFD_90")
  

  
  lifetables <- read_xlsx("lifetable/CBS probmort NL.xlsx", sheet = "input")
  lifetables$lograte <- log(lifetables$probmort)
  # Lifetables are rates so must be converted to probabilities
  lifetables$Overall <- 1 - exp(-lifetables$probmort)
  death_probability_nocomplications	<- data.frame(lifetables$age, lifetables$Overall)
  
  
  # Combine on log scale
  # Rows are for samples, columns are for ages
  death_probability_osteoporosis <- list()
  death_probability_osteoporosis$Overall <- 1 - exp(-exp(matrix(rep(input_parameters$death_log_hr_osteoporosis_mixedj, length(lifetables$lograte)) +   rep(lifetables$lograte, each = n_samples), nrow = n_samples)))
  colnames(death_probability_osteoporosis$Overall) <- lifetables$age
  
  
  
  
## INEFICIENT, but works, the below is to utilize parameters specific to children/adults whenever possible 
  
#from state to GIC
  #ages 2-18, note assumptions
    transition_matrices[, 1:17, "diagnosed no complications", "diagnosed GIC"] <- 
      transition_matrices[, 1:17, "diagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, 1:17, "diagnosed NHL", "diagnosed GIC"] <- 
      GIC_probability_GFD_all$GIC_probability_GFD_10
    transition_matrices[, 1:17, "Undiagnosed no complications", "Undiagnosed GIC"] <-
      transition_matrices[, 1:17, "Undiagnosed osteoporosis", "Undiagnosed GIC"] <-
      transition_matrices[, 1:17, "Undiagnosed NHL", "Undiagnosed GIC"] <-
      (1 - input_parameters$probability_late_diagnosis_children) * GIC_probability_noGFD_all$GIC_probability_noGFD_10
    transition_matrices[, 1:17, "Undiagnosed no complications", "diagnosed GIC"] <-
      transition_matrices[, 1:17, "Undiagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, 1:17, "Undiagnosed NHL", "diagnosed GIC"] <-
      input_parameters$probability_late_diagnosis_children * GIC_probability_noGFD_all$GIC_probability_noGFD_10
    transition_matrices[, 1:17,"Undiagnosed GIC", "diagnosed GIC"] <- input_parameters$probability_late_diagnosis_children 
    
    #ages 19-21, note assumptions
    transition_matrices[, 18:20, "diagnosed no complications", "diagnosed GIC"] <- 
      transition_matrices[, 18:20, "diagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, 18:20, "diagnosed NHL", "diagnosed GIC"] <- 
      GIC_probability_GFD_all$GIC_probability_GFD_20
    transition_matrices[, 18:20, "Undiagnosed no complications", "Undiagnosed GIC"] <-
      transition_matrices[, 18:20, "Undiagnosed osteoporosis", "Undiagnosed GIC"] <-
      transition_matrices[, 18:20, "Undiagnosed NHL", "Undiagnosed GIC"] <-
      (1 - input_parameters$probability_late_diagnosis_adults) * GIC_probability_noGFD_all$GIC_probability_noGFD_20
    transition_matrices[, 18:20, "Undiagnosed no complications", "diagnosed GIC"] <-
      transition_matrices[, 18:20, "Undiagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, 18:20, "Undiagnosed NHL", "diagnosed GIC"] <-
      input_parameters$probability_late_diagnosis_adults * GIC_probability_noGFD_all$GIC_probability_noGFD_20
    transition_matrices[, 18:20,"Undiagnosed GIC", "diagnosed GIC"] <- input_parameters$probability_late_diagnosis_adults   
    
  #ages 20-onwards
    n_agecategories <- (n_cycles/10) - 1
    for(i_age_category in c(2:n_agecategories)) {
    transition_matrices[, (c(1:10) + i_age_category * 10), "diagnosed no complications", "diagnosed GIC"] <- 
      transition_matrices[, (c(1:10) + i_age_category * 10), "diagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, (c(1:10) + i_age_category * 10), "diagnosed NHL", "diagnosed GIC"] <- 
      GIC_probability_GFD_all[, starting_age_column + i_age_category]
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "Undiagnosed GIC"] <-
      transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed osteoporosis", "Undiagnosed GIC"] <-
      transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed NHL", "Undiagnosed GIC"] <-
      (1 - input_parameters$probability_late_diagnosis_adults) * GIC_probability_noGFD_all[, starting_age_column + i_age_category] 
    transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "diagnosed GIC"] <-
      transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed osteoporosis", "diagnosed GIC"] <-
      transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed NHL", "diagnosed GIC"] <-
      input_parameters$probability_late_diagnosis_adults * GIC_probability_noGFD_all[, starting_age_column + i_age_category] 
    transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed GIC", "diagnosed GIC"] <- input_parameters$probability_late_diagnosis_adults 
    }
  
    
    
    
    #from state to NHL
    #ages 2-18, note assumptions
    transition_matrices[, 1:17, "diagnosed no complications", "diagnosed NHL"] <- 
      transition_matrices[, 1:17, "diagnosed osteoporosis", "diagnosed NHL"] <-
      NHL_probability_GFD_all$NHL_probability_GFD_10
    transition_matrices[, 1:17, "Undiagnosed no complications", "Undiagnosed NHL"] <-
      transition_matrices[, 1:17, "Undiagnosed osteoporosis", "Undiagnosed NHL"] <-
      (1 - input_parameters$probability_late_diagnosis_children) * NHL_probability_noGFD_all$NHL_probability_noGFD_10
    transition_matrices[, 1:17, "Undiagnosed no complications", "diagnosed NHL"] <-
      transition_matrices[, 1:17, "Undiagnosed osteoporosis", "diagnosed NHL"] <-
      input_parameters$probability_late_diagnosis_children * NHL_probability_noGFD_all$NHL_probability_noGFD_10
    transition_matrices[, 1:17,"Undiagnosed NHL", "diagnosed NHL"] <- input_parameters$probability_late_diagnosis_children - 
      (transition_matrices[, 1:17, "Undiagnosed NHL","diagnosed GIC"])
    
    #ages 19-21, note assumptions
    transition_matrices[, 18:20, "diagnosed no complications", "diagnosed NHL"] <- 
      transition_matrices[, 18:20, "diagnosed osteoporosis", "diagnosed NHL"] <-
      NHL_probability_GFD_all$NHL_probability_GFD_20
    transition_matrices[, 18:20, "Undiagnosed no complications", "Undiagnosed NHL"] <-
      transition_matrices[, 18:20, "Undiagnosed osteoporosis", "Undiagnosed NHL"] <-
      (1 - input_parameters$probability_late_diagnosis_adults) * NHL_probability_noGFD_all$NHL_probability_noGFD_20
    transition_matrices[, 18:20, "Undiagnosed no complications", "diagnosed NHL"] <-
      transition_matrices[, 18:20, "Undiagnosed osteoporosis", "diagnosed NHL"] <-
      input_parameters$probability_late_diagnosis_adults * NHL_probability_noGFD_all$NHL_probability_noGFD_20
    transition_matrices[, 18:20,"Undiagnosed NHL", "diagnosed NHL"] <- input_parameters$probability_late_diagnosis_adults - 
      (transition_matrices[, 18:20, "Undiagnosed NHL","diagnosed GIC"])   
    
    #ages 20-onwards
    n_agecategories <- (n_cycles/10) - 1
    for(i_age_category in c(2:n_agecategories)) {
      transition_matrices[, (c(1:10) + i_age_category * 10), "diagnosed no complications", "diagnosed NHL"] <- 
        transition_matrices[, (c(1:10) + i_age_category * 10), "diagnosed osteoporosis", "diagnosed NHL"] <-
        NHL_probability_GFD_all[, starting_age_column + i_age_category]
      transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "Undiagnosed NHL"] <-
        transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed osteoporosis", "Undiagnosed NHL"] <-
        (1 - input_parameters$probability_late_diagnosis_adults) * NHL_probability_noGFD_all[, starting_age_column + i_age_category] 
      transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "diagnosed NHL"] <-
        transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed osteoporosis", "diagnosed NHL"] <-
        input_parameters$probability_late_diagnosis_adults * NHL_probability_noGFD_all[, starting_age_column + i_age_category] 
      transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed NHL", "diagnosed NHL"] <- input_parameters$probability_late_diagnosis_adults - 
        (transition_matrices[,(c(1:10) + i_age_category * 10),"Undiagnosed NHL","diagnosed GIC"])
    }
    
    
    
    
    
  
    
  # from state to OSTEOPOROSIS  
    #ages 2-18
    transition_matrices[, 1:17, "diagnosed no complications", "diagnosed osteoporosis"] <- osteoporosis_probability_GFD_all$osteoporosis_probability_GFD_10
    transition_matrices[, 1:17, "Undiagnosed no complications", "Undiagnosed osteoporosis"] <- (1 - input_parameters$probability_late_diagnosis_children) * osteoporosis_probability_noGFD_all$osteoporosis_probability_noGFD_10
    transition_matrices[, 1:17, "Undiagnosed no complications", "diagnosed osteoporosis"] <-  input_parameters$probability_late_diagnosis_children * osteoporosis_probability_noGFD_all$osteoporosis_probability_noGFD_10
    transition_matrices[, 1:17,"Undiagnosed no complications", "diagnosed no complications"] <-
      input_parameters$probability_late_diagnosis_children -
      (transition_matrices[, 1:17, "Undiagnosed no complications", "diagnosed osteoporosis"] +
         (transition_matrices[, 1:17, "Undiagnosed no complications", "diagnosed NHL"] ) +
         (transition_matrices[, 1:17, "Undiagnosed no complications", "diagnosed GIC"] ) ) 
    transition_matrices[,1:17 ,"Undiagnosed osteoporosis", "diagnosed osteoporosis"] <-
      input_parameters$probability_late_diagnosis_children -
      ((transition_matrices[, 1:17,"Undiagnosed osteoporosis", "diagnosed NHL"]) +
         (transition_matrices[, 1:17,"Undiagnosed osteoporosis", "diagnosed GIC"]))
    
    #ages 19-21
    
    transition_matrices[, 18:20, "diagnosed no complications", "diagnosed osteoporosis"] <- osteoporosis_probability_GFD_all$osteoporosis_probability_GFD_20
    transition_matrices[, 18:20, "Undiagnosed no complications", "Undiagnosed osteoporosis"] <- (1 - input_parameters$probability_late_diagnosis_adults) * osteoporosis_probability_noGFD_all$osteoporosis_probability_noGFD_20
    transition_matrices[, 18:20, "Undiagnosed no complications", "diagnosed osteoporosis"] <-  input_parameters$probability_late_diagnosis_adults * osteoporosis_probability_noGFD_all$osteoporosis_probability_noGFD_20
    transition_matrices[, 18:20,"Undiagnosed no complications", "diagnosed no complications"] <-
      input_parameters$probability_late_diagnosis_adults -
      (transition_matrices[, 18:20, "Undiagnosed no complications", "diagnosed osteoporosis"] +
         (transition_matrices[, 18:20, "Undiagnosed no complications", "diagnosed NHL"] ) +
         (transition_matrices[, 18:20, "Undiagnosed no complications", "diagnosed GIC"] ) ) 
    transition_matrices[,18:20 ,"Undiagnosed osteoporosis", "diagnosed osteoporosis"] <-
      input_parameters$probability_late_diagnosis_adults -
      ((transition_matrices[, 18:20,"Undiagnosed osteoporosis", "diagnosed NHL"]) +
         (transition_matrices[, 18:20,"Undiagnosed osteoporosis", "diagnosed GIC"]))
  
    #ages 20-onwards
   n_agecategories <- (n_cycles/10) - 1
   for(i_age_category in c(2:n_agecategories)) {
     transition_matrices[, (c(1:10) + i_age_category * 10), "diagnosed no complications", "diagnosed osteoporosis"] <- osteoporosis_probability_GFD_all[, starting_age_column + i_age_category]
     transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "Undiagnosed osteoporosis"] <- (1 - input_parameters$probability_late_diagnosis_adults) * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
     transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "diagnosed osteoporosis"] <-  input_parameters$probability_late_diagnosis_adults * osteoporosis_probability_noGFD_all[, starting_age_column + i_age_category]
     transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed no complications", "diagnosed no complications"] <-
          input_parameters$probability_late_diagnosis_adults -
          (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "diagnosed osteoporosis"] +
          (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "diagnosed NHL"] ) +
          (transition_matrices[, (c(1:10) + i_age_category * 10), "Undiagnosed no complications", "diagnosed GIC"] ) ) 
     transition_matrices[,(c(1:10) + i_age_category * 10) ,"Undiagnosed osteoporosis", "diagnosed osteoporosis"] <-
          input_parameters$probability_late_diagnosis_adults -
          ((transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed osteoporosis", "diagnosed NHL"]) +
          (transition_matrices[, (c(1:10) + i_age_category * 10),"Undiagnosed osteoporosis", "diagnosed GIC"]))
   }
   
                                       
  n_ages <- 82 - starting_age
  
  # NHL and osteoporosis death probabilities very high so instead of assuming no competing risks we assume
  # it reduces all other transitions (e.g., the 33% who do not die of NHL are divided among remaining probabilities)
  # This avoids negative probabilities of remaining in the same state
  death_state_index <- which(state_names == "Death")
  for (i_age in 1:n_ages){
    transition_matrices[, i_age, "diagnosed no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "diagnosed osteoporosis", "Death"] <- death_probability_osteoporosis$Overall[, starting_age + i_age]
    transition_matrices[, i_age, "diagnosed osteoporosis", -death_state_index] <-  transition_matrices[, i_age, "diagnosed osteoporosis", -death_state_index] * (1 - death_probability_osteoporosis$Overall[, starting_age + i_age])
    transition_matrices[, i_age, "Undiagnosed no complications", "Death"] <- death_probability_nocomplications[starting_age + i_age, 2]
    transition_matrices[, i_age, "Undiagnosed osteoporosis", "Death"] <- death_probability_osteoporosis$Overall[, starting_age + i_age]
    transition_matrices[, i_age, "Undiagnosed osteoporosis", -death_state_index] <- transition_matrices[, i_age, "Undiagnosed osteoporosis", -death_state_index] * (1 - death_probability_osteoporosis$Overall[, starting_age + i_age])
  }
  
  # Probabilities of death form NHL
  transition_matrices[, , "diagnosed NHL", "Death"] <- 
    transition_matrices[, , "Undiagnosed NHL", "Death"] <- input_parameters$death_probability_NHL
  
  # Probabilities of death form GIC
  transition_matrices[, , "diagnosed GIC", "Death"] <- 
    transition_matrices[, , "Undiagnosed GIC", "Death"] <- input_parameters$death_probability_GIC
  
  # Multiply other transitions to prevent probabilities exceeding 1  
  for(i_sample in 1:n_samples) {
    transition_matrices[i_sample, , "Undiagnosed NHL", -death_state_index] <- 
      transition_matrices[i_sample, , "Undiagnosed NHL", -death_state_index] * (1 - input_parameters$death_probability_NHL[i_sample])
    transition_matrices[i_sample, , "diagnosed NHL", -death_state_index] <- 
      transition_matrices[i_sample, , "diagnosed NHL", -death_state_index] * (1 - input_parameters$death_probability_NHL[i_sample])
  }
  
  for(i_sample in 1:n_samples) {
    transition_matrices[i_sample, , "Undiagnosed GIC", -death_state_index] <- 
      transition_matrices[i_sample, , "Undiagnosed GIC", -death_state_index] * (1 - input_parameters$death_probability_GIC[i_sample])
    transition_matrices[i_sample, , "diagnosed GIC", -death_state_index] <- 
      transition_matrices[i_sample, , "diagnosed GIC", -death_state_index] * (1 - input_parameters$death_probability_GIC[i_sample])
  }
  
  
  #Complete the matrix by adding the complement of all probabilities  
  for(i_state in 1:length(state_names)) {
    transition_matrices[, , i_state, i_state] <- 1 - 
      apply(transition_matrices[, , i_state, -i_state], c(1,2), sum, na.rm=TRUE)
  }
  
  
  transition_matrices[, , , ] [is.na(transition_matrices[, , , ] )] <- 0
  
  # Check that rows sum to 1
   rowSums (transition_matrices[1, 4, , ], na.rm = FALSE , dims = 1)
  
  return(transition_matrices[, , , ])
}

