
#project: Glutenscreen markov
#action: main model
#title: glutenscreen_main.R
#version: g

# Load BCEA library to help analyse and visualise results
pacman::p_load(
  BCEA,
  readxl,
  ggplot2,
  reshape,
  dplyr,
  tictoc,
  SimDesign,
  BCEA,
  dplyr,
  EnvStats,
  Rmisc,
  car,
  Rcpp)

set.seed(14143)



source('generate_input_parameters_e.R')
source('generate_transition_matrices_b.R')
source('generate_state_qalys.R')
source('generate_state_costs_b.R')
source('convert_transition_matrices_to_df.R')
source('convert_cohort_vectors_to_df.R')
source('generate_net_benefit_c.R')


test <- "POCT"


# Define global simulation parameters
n_samples <- 1000
                                                                                      
n_states <- 9
state_names <- c("diagnosed no complications",
                 "diagnosed osteoporosis",
                 "diagnosed NHL",
                 "diagnosed GIC",
                 "Undiagnosed no complications",
                 "Undiagnosed osteoporosis",
                 "Undiagnosed NHL",
                 "Undiagnosed GIC",
                 "Death")


starting_age <- 2 

n_cycles <- 82-starting_age 

lambda <- 20000 #WTP threshold

#define strategies (excel spreadsheet 'Model info and specifications' sheet 'Quest_test_sens_spec' displays the same process)
constant <- 0.01  #this represents a prevalence. As long as this nr is not 1 or 0 (implying 100 or 0% prevalence), then we get the same values for sens_riskfactor and spec_riskfactor in the end (see below) 
response_rate <- c(0.8) #assumed response rate to the questionnaire
sens_q <- c(0.9, 0.8, 0.7) #realistic sensitivities for the questionnaire 
spec_q <- c(0.3, 0.5, 0.7) #realistic specificities for the questionnaire
test_rate <- c(0.9) #of those that are eligible for testing with POCT after answering the questionnaire (i.e. true positives + false positives based on the above sens_q and spec_q), a small portion is lost (e.g. parental consent, no show, etc)

combinations <- expand.grid(constant = constant,
                            response_rate = response_rate,
                            sens_q = sens_q,
                            spec_q = spec_q,
                            test_rate = test_rate)

#calculate 'effective' sens/spec (synonymous to Bristol's sens_riskfactor and spec_riskfactor) taking into account the response rate, sens/spec of questionnaire, and test rate
combinations$true_cd <- combinations$constant
combinations$true_no_cd <- 1-combinations$constant
combinations$cd_responded <- combinations$true_cd*combinations$response_rate
combinations$cd_not_responded <- combinations$true_cd- combinations$cd_responded
combinations$no_cd_responded <- combinations$true_no_cd* combinations$response_rate
combinations$no_cd_not_responded <- combinations$true_no_cd - combinations$no_cd_responded

combinations$tp_q <- combinations$cd_responded* combinations$sens_q
combinations$fn_q <- combinations$cd_responded- combinations$tp_q
combinations$tn_q <- combinations$no_cd_responded* combinations$spec_q
combinations$fp_q <- combinations$no_cd_responded - combinations$tn_q

combinations$cd_tested <- combinations$tp_q * combinations$test_rate
combinations$cd_not_tested <- combinations$fn_q + ( combinations$tp_q - combinations$cd_tested)
combinations$no_cd_tested <- combinations$fp_q * combinations$test_rate
combinations$no_cd_not_tested <-  combinations$tn_q + ( combinations$fp_q - combinations$no_cd_tested)

combinations$tp_riskfactor <- combinations$cd_tested #necessary to calculate the final sens_riskfactor and spec_riskfactor values to be used in the generate_net_benefit function
combinations$fn_riskfactor <- combinations$true_cd - combinations$tp_riskfactor #necessary to calculate the final sens_riskfactor and spec_riskfactor values to be used in the generate_net_benefit function
combinations$fp_riskfactor <- combinations$no_cd_tested #necessary to calculate the final sens_riskfactor and spec_riskfactor values to be used in the generate_net_benefit function
combinations$tn_riskfactor <- combinations$no_cd_not_tested + combinations$no_cd_not_responded #necessary to calculate the final sens_riskfactor and spec_riskfactor values to be used in the generate_net_benefit function

combinations$sens_riskfactor <- combinations$tp_riskfactor/(combinations$tp_riskfactor + combinations$fn_riskfactor)
combinations$spec_riskfactor <- combinations$tn_riskfactor/(combinations$tn_riskfactor+ combinations$fp_riskfactor)
combinations$prop_cd_missed <- combinations$fn_riskfactor/combinations$constant
combinations$ppv <- combinations$tp_riskfactor/ (combinations$tp_riskfactor+ combinations$fp_riskfactor)

comb_ref <- combinations #just for future reference
combinations <- combinations[,24:25]

#add Bristol combinations of interest, for potential relevance to our setting
brismod_interest <- as.data.frame(read.csv("brismod/data/included_strategies.csv"))
brismod_interest <- brismod_interest[1:6, 8:9]
colnames(brismod_interest) <- c("sens_riskfactor", "spec_riskfactor")
combinations <- rbind(combinations, brismod_interest)
combinations$sens_riskfactor[combinations$sens_riskfactor == 1] <- 0.9999


# Codes and numbers of all risk prediction strategies 
combinations$x <- paste(combinations$sens_riskfactor, combinations$spec_riskfactor)
combinations_names <- combinations$x
n_combinations <- length(combinations$x)


# Define the number and names of tests

n_tests <- (1 * n_combinations) + 1

t_names <-  c("No screening", outer(combinations_names, test, FUN = "paste")[1:n_tests-1])




# Generate the input parameters
# This will be converted into transition matrix, state costs, and state utilities
input_parameters <- generate_input_parameters_e(n_samples = n_samples,
                                                starting_age = starting_age)


transition_matrices <- generate_transition_matrices(input_parameters)

# Run the Markov model to get the model outputs
output <- generate_net_benefit_c(input_parameters = input_parameters, 
                                 strategies_of_interest = strategies_of_interest, 
                                 transition_matrices = transition_matrices,
                                 combinations = combinations,
                                 lambda=20000)

holder <- c("No screening", "0.576 0.496 POCT", "0.648 0.784 POCT", "0.9999 0 POCT")
#holder <- t_names[! t_names %in% c("0.9999 0 POCT" )] 
         # "0.882 0.417 POCT",
         # "0.807 0.61 POCT"  ,
         # "0.667 0.872 POCT"  ,
         # "0.533 0.952 POCT"  ,
         # "0.331 0.987 POCT"   )]
              
                            
m <- bcea(e = t(output$total_qalys[holder,]), 
          c = t(output$total_costs[holder,]), ref = 1, 
          interventions = holder)
summary(m)

ceac.plot(m, graph = c("base"), pos = "topright",
          cex.main = 0.2,
          line = list(colors = c(1:length(holder)))) #all two-way comparisons

#multiple comparisons
mce <- multi.ce(m)

#pdf(file = paste0("results/_ceac_all.pdf"))
ceac.plot(mce, graph = c("base"), pos = "topright",
          cex.main = 0.2,
          line = list(colors = c(1:length(holder)))) #CEAC 
#dev.off()

############################################################################################                            

strategies_excluded <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit < 0)) #strategies with ENB less than no screening
strategies_included <- names(subset(output$incremental_net_benefit,output$incremental_net_benefit > 0)) #strategies with ENB greater than no screening




write.csv(t(output$total_qalys), paste0("results/total qalys.csv"))
write.csv(t(output$total_costs), paste0("results/total costs.csv"))

#table of costs and QALYs
format_results<-function(x, n_digits = 2) {
  paste(round(mean(x, na.rm = TRUE),digits = n_digits), " (",
        round(quantile(x, probs = 0.025, na.rm = TRUE), digits = n_digits), ", ", 
        round(quantile(x, probs = 0.975, na.rm = TRUE), digits = n_digits),")",sep="")
}

################################################################################################################
# This where specific strategies are selected ##JH CHANGES: from 'strategies_of_interest' to 't_names'###########

t_names_qalys <- output$total_qalys[t_names, ]
t_names_costs <- output$total_costs[t_names, ]
t_names_inetbenefit <- output$all_incremental_net_benefit[t_names, ]
t_names_qalys_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_costs_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_inetbenefit_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))

for (i in 1:length(t_names)) { 
  t_names_qalys_table[, i] <- format_results(t_names_qalys[i, ], n_digits = 2)
  t_names_costs_table[, i] <- format_results(t_names_costs[i, ], n_digits = 0)
  t_names_inetbenefit_table[, i] <- format_results(t_names_inetbenefit[i, ], n_digits = 0)
}

# Create a table comparing costs, QALYs, and incrmental net benefit for strategies of interest
results_table <- data.frame(
  c("", rep(combinations$sens_riskfactor, length(n_tests))),
  c("", rep(combinations$spec_riskfactor, length(n_tests))),
  t(t_names_costs_table), t(t_names_qalys_table), t(t_names_inetbenefit_table))
colnames(results_table) <- c("Sensitivity", "Specificity", "Costs", "QALYs", "Incremental net benefit v no screening")

# Change the sensitivity and specificity of 0.9999 back to 1 for ease of presentation
results_table$Sensitivity[results_table$Sensitivity == 0.9999] <- 1
results_table$Specificity[results_table$Specificity == 0.9999] <- 1
rownames(results_table) <- t_names
write.csv(results_table, paste0("results/table of results.csv"))

#cost breakdown
t_names_testcosts <- output$test_costs_applied[, t_names]
t_names_diagnosiscosts <- output$diagnosis_costs[, t_names]

t_names_cyclecosts <- output$cycle_costs[t_names]
t_names_testcosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_diagnosiscosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_cyclecosts_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
for (i in 1:length(t_names)) { 
  t_names_testcosts_table[, i] <- format_results(t_names_testcosts[, i])
  t_names_diagnosiscosts_table[, i] <- format_results(t_names_diagnosiscosts[, i])
  t_names_cyclecosts_table[i] <- format_results(t_names_cyclecosts[i])
}

cost_breakdown_table <- data.frame(t(t_names_testcosts_table), t(t_names_diagnosiscosts_table), t(t_names_cyclecosts_table))
colnames(cost_breakdown_table) <- c("Test costs", "Diagnosis costs", "Cycle costs")
rownames(cost_breakdown_table) <- t_names

write.csv(cost_breakdown_table, paste0("results/cost breakdown table.csv"))


#utility breakdown
t_names_disutilitybiopsy <- output$disutility_biopsy[, t_names]
t_names_disutilitybiopsywait <- output$disutility_biopsy_wait[, t_names]

t_names_cycleqalys <- output$cycle_qalys[t_names]
t_names_disutilitybiopsy_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_disutilitybiopsywait_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
t_names_cycleqalys_table <- array(dim=c(1, length(t_names)), dimnames = list(NULL, NULL))
for (i in 1:length(t_names)) { 
  t_names_disutilitybiopsy_table[, i] <- format_results(t_names_disutilitybiopsy[, i], n_digits = 4)
  t_names_disutilitybiopsywait_table[, i] <- format_results(t_names_disutilitybiopsywait[, i], n_digits = 4)
  t_names_cycleqalys_table[i] <- format_results(t_names_cycleqalys[i], n_digits = 4)
}

qaly_breakdown_table <- data.frame(t(t_names_cycleqalys_table), t(t_names_disutilitybiopsy_table), t(t_names_disutilitybiopsywait_table))
colnames(qaly_breakdown_table) <- c("Cycle QALYs", "Disutility biopsy", "Disutility biopsy wait")
rownames(qaly_breakdown_table) <- t_names

write.csv(qaly_breakdown_table, paste0("results/qaly breakdown table.csv"))

#time in states
time_in_states <- output$time_in_states

write.csv(time_in_states, paste0("results/time in states.csv"))




