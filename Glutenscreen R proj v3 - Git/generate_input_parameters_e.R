
#project: Glutenscreen markov
#action: generate input parameters function
#title: generate_input_parameters.R


# Utility functions
expit <- function(logO) {
  return(exp(logO)/(1 + exp(logO)))
}
logit <- function(p) {
  return(log(p/(1-p)))
}

generate_input_parameters_e <- function(n_samples,
                                      starting_age = starting_age) {
  
  
  
  duration_of_symptoms_adults <- 10.93   #Violato et al 
  duration_of_symptoms_adults_sd <- 13.10 
  duration_of_symptoms_adults_location <- log(duration_of_symptoms_adults ^ 2 / sqrt(duration_of_symptoms_adults_sd ^ 2 + duration_of_symptoms_adults ^ 2))
  duration_of_symptoms_adults_shape <- sqrt(log(1 + (duration_of_symptoms_adults_sd ^ 2 / duration_of_symptoms_adults ^ 2)))
  duration_of_symptoms_adults <- rlnorm(n = n_samples, duration_of_symptoms_adults_location,  duration_of_symptoms_adults_shape)     #calculated from Violato et al 2019
  rate_of_symptoms_adults <- 1 / duration_of_symptoms_adults
  probability_late_diagnosis_adults <- 1 - exp(-rate_of_symptoms_adults)

  
  duration_of_symptoms_children <-  3.34  #Violato et al
  duration_of_symptoms_children_sd <- 3.71 
  duration_of_symptoms_children_location <- log(duration_of_symptoms_children ^ 2 / sqrt(duration_of_symptoms_children_sd ^ 2 + duration_of_symptoms_children ^ 2))
  duration_of_symptoms_children_shape <- sqrt(log(1 + (duration_of_symptoms_children_sd ^ 2 / duration_of_symptoms_children ^ 2)))
  duration_of_symptoms_children <- rlnorm(n = n_samples, duration_of_symptoms_children_location,  duration_of_symptoms_children_shape)     #calculated from Violato et al 2019
  rate_of_symptoms_children <- 1 / duration_of_symptoms_children
  probability_late_diagnosis_children <- 1 - exp(-rate_of_symptoms_children)
  
  # duration_of_symptoms_children_2 <- rtri(n=n_samples,
  #                                       mode=(87/12),
  #                                       min=0.0001,
  #                                       max=(637.7/12)) #perhaps more realistic. From https://www.sciencedirect.com/science/article/pii/S1590865816304753 ...
  # rate_of_symptoms_children_2 <- 1 / duration_of_symptoms_children_2
  # probability_late_diagnosis_children_2 <- 1 - exp(-rate_of_symptoms_children_2)

  # Use the appropriate prevalence
  prevalence = as.data.frame(readxl::read_excel(path = "prevalence/adapted from CPRD prevalence.xlsx", sheet = "mixed")) 
  
  #Initial cohort at diagnosis - depends on age at diagnosis 
  # e.g. 10 year old starts with 10-20 prevalence and 18 year old starts with 10-20 prevalence. 
  starting_age_category <- 1
  # prevalence$Age.categories == starting_age
  probability_osteoporosis <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "Osteoporosis_r"], shape2 = prevalence[starting_age_category, "N"] - prevalence[starting_age_category, "Osteoporosis_r"])
  probability_NHL <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "NHL_r"], shape2 = prevalence[starting_age_category, "N"] - prevalence[starting_age_category, "NHL_r"])
  probability_GIC <- rbeta(n=n_samples, shape1 = prevalence[starting_age_category, "GIC_r"], shape2 = prevalence[starting_age_category, "GIC_N"] - prevalence[starting_age_category, "GIC_r"])
  probability_nocomplications <- 1 - probability_osteoporosis - probability_NHL - probability_GIC
  
  # IDA prevalence changes with age of cohort so is age stratified
  # Prevalence of osteoporosis and NHL are combined with incidence rates to model prevalence changing with age
  probability_IDA_0 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[1], shape2 = (prevalence$N[1] - prevalence$IDA_r[1]))
  probability_IDA_10 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[2], shape2 = (prevalence$N[2] - prevalence$IDA_r[2]))
  probability_IDA_20 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[3], shape2 = (prevalence$N[3] - prevalence$IDA_r[3]))
  probability_IDA_30 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[4], shape2 = (prevalence$N[4] - prevalence$IDA_r[4]))
  probability_IDA_40 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[5], shape2 = (prevalence$N[5] - prevalence$IDA_r[5]))
  probability_IDA_50 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[6], shape2 = (prevalence$N[6] - prevalence$IDA_r[6]))
  probability_IDA_60 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[7], shape2 = (prevalence$N[7] - prevalence$IDA_r[7]))
  probability_IDA_70 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[8], shape2 = (prevalence$N[8] - prevalence$IDA_r[8]))
  probability_IDA_80 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[9], shape2 = (prevalence$N[9] - prevalence$IDA_r[9]))
  probability_IDA_90 <- rbeta(n=n_samples, shape1 = prevalence$IDA_r[10], shape2 = (prevalence$N[10] - prevalence$IDA_r[10]))
  probability_IDA <- data.frame(probability_IDA_0 , probability_IDA_10 , probability_IDA_20 , probability_IDA_30 , probability_IDA_40
                                ,probability_IDA_50 , probability_IDA_60 , probability_IDA_70 , probability_IDA_80 , probability_IDA_90)
  
  
  #####Developing consequence
  
  ## Developing Osteoporosis (Osteoporosis probabilities On GFD)
  
  Osteoporosis_probability <- as.data.frame(readxl::read_excel(path = "baseline_rates/ostp_rates_NL.xlsx", sheet = "Sheet1")) 
  Osteoporosis_probability$Osteoporosis_rate <- Osteoporosis_probability$Osteoporosis_rate
  Osteoporosis_probability$Osteoporosis_rate_low <- 0.75*Osteoporosis_probability$Osteoporosis_rate
  Osteoporosis_probability$Osteoporosis_rate_high <- 1.25*Osteoporosis_probability$Osteoporosis_rate

  
  
# Log odds ratio for diagnosed CD
  log_or_osteoporosis_GFD <- rnorm(n_samples, mean = log(1.43), sd = ((log(1.78) - log(1.15))/(2*1.96)))  #Olmos 2008
  log_or_osteoporosis_noGFD <- rnorm(n_samples, mean = log(1.77), sd = ((log(2.65) - log(1.18))/(2*1.96))) #Hujoel 2018
  
  

  # Log rates in general population
  osteoporosis_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(osteoporosis_lograte) <- paste0("osteoporosis_lograte_", seq(0, 90, 10))

  for(i_age_category in 1:10) {
    osteoporosis_lograte[, i_age_category] <- rtri(n=n_samples,
                                                   mode = log(Osteoporosis_probability$Osteoporosis_rate[i_age_category]),
                                                   min =  log(Osteoporosis_probability$Osteoporosis_rate_low[i_age_category]),
                                                   max =  log(Osteoporosis_probability$Osteoporosis_rate_high[i_age_category])
                                                   )
  }
  
  ##developing NHL
  
  NHL_probability <- as.data.frame(readxl::read_excel(path = "baseline_rates/NHL_rates.xlsx", sheet = "input"))
  NHL_probability$NHL_rate_low <- 0.75*NHL_probability$NHL_rate
  NHL_probability$NHL_rate_high <- 1.25*NHL_probability$NHL_rate
  
  # Log incidence ratios for NHL
  log_rr_NHL_GFD <- rnorm(n_samples, mean = log(3.28), sd = (log(6.28) - log(1.49))/(2*1.96))
  log_rr_NHL_noGFD <- rnorm(n_samples, mean = log(4.7), sd = (log(7.3) - log(2.9))/(2*1.96))
  
  # Log rates in general population
  NHL_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(NHL_lograte) <- paste0("NHL_lograte_", seq(0, 90, 10))
  for(i_age_category in 1:10) {
    NHL_lograte[, i_age_category] <- rtri(n=n_samples, 
                                          mode = log(NHL_probability$NHL_rate[i_age_category]),
                                          min = log(NHL_probability$NHL_rate_low[i_age_category]) ,
                                          max = log(NHL_probability$NHL_rate_high[i_age_category])
                                          )
  }
  
  
  ##developing GIC 
  
  GIC_probability <- as.data.frame(readxl::read_excel(path = "baseline_rates/GIC_rates.xlsx", sheet = "input"))
  GIC_probability$GIC_rate_low <- 0.75*GIC_probability$GIC_rate
  GIC_probability$GIC_rate_high <- 1.25*GIC_probability$GIC_rate
  #Ensure no rates are zero
  
  
  # Log incidence ratios for GIC from studies
   log_rr_GIC_GFD <- rnorm(n_samples, mean = log(1.32), sd = (log(1.43) - log(1.22))/(2*1.96))
   log_rr_GIC_noGFD <- rnorm(n_samples, mean = log(2.33), sd = (log(4.04) - log(1.35))/(2*1.96))
  
  
  
  # Log rates in general population
  GIC_lograte <- matrix(NA, nrow = n_samples, ncol = 10)
  colnames(GIC_lograte) <- paste0("GIC_lograte_", seq(0, 90, 10))
  for(i_age_category in 1:10) {
    GIC_lograte[, i_age_category] <- rtri(n=n_samples, 
                                          mode = log(GIC_probability$GIC_rate[i_age_category]),
                                          min =  log(GIC_probability$GIC_rate_low[i_age_category]),
                                          max = log(GIC_probability$GIC_rate_high[i_age_category])
                                          )
  }
  
  
  ##Dying from consequence
  death_log_hazard_NHL <- exp(-2.6026)  # 2.6026 based on IKNL (indolent and aggresive)
  death_probability_NHL <-	rtri( n=n_samples,
                                  min= 0.8*(1-exp(-exp(death_log_hazard_NHL))),
                                  max= 1.2*(1-exp(-exp(death_log_hazard_NHL))),
                                  mode = 1-exp(-exp(death_log_hazard_NHL))   
                                  )

  
  # Death OSTP - from Klop 2014 (other sources: Farahmand 2005, Abrahamsen 2009)
  death_log_hr_osteoporosis_mixedj <- rnorm(n = n_samples, mean = log(2.80),
                                            sd = (log(3.04) - log(2.58))/(2*1.96))
  
  # Death GIC - from IKNL (spreadsheet "GIC rates death"), using method in Bristol
  death_log_hazard_GIC <- exp(-0.9729)
  death_probability_GIC <- rtri(n=n_samples,
                                mode = 1-exp(-exp(death_log_hazard_GIC)), 
                                min = 0.8*(1-exp(-exp(death_log_hazard_GIC))),
                                max = 1.2*(1-exp(-exp(death_log_hazard_GIC)))
  )
 
  
  #############################################################################
  ## Utilities ################################################################
  #############################################################################
  ##CHECK FOR UPDATES
  
  utility_GFD_adults <- 0.85  
  utility_GFdse_adults <- ((0.86-0.84)/3.92) 
  utility_GFDalpha_adults <- (utility_GFD_adults ^ 2 * (1 - utility_GFD_adults)/utility_GFdse_adults ^ 2) - utility_GFD_adults
  utility_GFDbeta_adults <- (utility_GFDalpha_adults / utility_GFD_adults) - utility_GFDalpha_adults
  utility_GFD_adults <- rbeta(n = n_samples, shape1 = utility_GFDalpha_adults, shape2 = utility_GFDbeta_adults)
  
  utility_GFD_children <- 0.88
  utility_GFdse_children <- ((0.92-0.85)/3.92)
  utility_GFDalpha_children <- (utility_GFD_children ^ 2 * (1 - utility_GFD_children)/utility_GFdse_children ^ 2) - utility_GFD_children
  utility_GFDbeta_children <- (utility_GFDalpha_children / utility_GFD_children) - utility_GFDalpha_children
  utility_GFD_children <- rbeta(n = n_samples, shape1 = utility_GFDalpha_children, shape2 = utility_GFDbeta_children)
  
  utility_undiagnosedCD_adults <-  0.65
  utility_undiagnosedCD_se_adults <- (0.67 - 0.63)/3.92
  utility_undiagnosedCD_alpha_adults <- (utility_undiagnosedCD_adults ^ 2 * (1 - utility_undiagnosedCD_adults)/utility_undiagnosedCD_se_adults ^ 2) - utility_undiagnosedCD_adults
  utility_undiagnosedCD_beta_adults <- (utility_undiagnosedCD_alpha_adults/utility_undiagnosedCD_adults) - utility_undiagnosedCD_alpha_adults
  utility_undiagnosedCD_adults <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha_adults, shape2 = utility_undiagnosedCD_beta_adults)

  utility_undiagnosedCD_children <-  0.65
  utility_undiagnosedCD_se_children <- (0.67 - 0.63)/3.92
  utility_undiagnosedCD_alpha_children <- (utility_undiagnosedCD_children ^ 2 * (1 - utility_undiagnosedCD_children)/utility_undiagnosedCD_se_children ^ 2) - utility_undiagnosedCD_children
  utility_undiagnosedCD_beta_children <- (utility_undiagnosedCD_alpha_children/utility_undiagnosedCD_children) - utility_undiagnosedCD_alpha_children
  utility_undiagnosedCD_children <- rbeta(n = n_samples, shape1 = utility_undiagnosedCD_alpha_children, shape2 = utility_undiagnosedCD_beta_children)

  
  #####
  
  probability_hipfracture <- 151/100000 # in NL, from Lotters et al 2016
  probability_vertebralfracture <- 53/100000 # in NL, from Lotters et al 2016
  probability_wristfracture <- 146/100000 # in NL, from Lotters et al 2016
  disutility_hipfracture <-   0.839 - 0.59  #0.839 corresponds to Dutch eq5d5l norm from Versteegh (SD: 0.183) ages 60-70
  disutility_hipfractureSE <- (((  0.839 - (1.1*0.59)) - (  0.839 - (0.9*0.59)))/3.92)
  disutility_hipfracture_alpha <- (disutility_hipfracture ^ 2 * (1 - disutility_hipfracture)/disutility_hipfractureSE ^ 2) - disutility_hipfracture
  disutility_hipfracture_beta <- (disutility_hipfracture_alpha/disutility_hipfracture) - disutility_hipfracture_alpha
  disutility_hipfracture <- rbeta(n = n_samples, shape1 = disutility_hipfracture_alpha, shape2 = disutility_hipfracture_beta)
  disutility_wristfracture <-   0.839 - 0.55
  disutility_wristfractureSE <- (((  0.839 - (1.1*0.55)) - (  0.839 - (0.9*0.55)))/3.92)
  disutility_wristfracture_alpha <- (disutility_wristfracture ^ 2 * (1 - disutility_wristfracture)/disutility_wristfractureSE ^ 2) - disutility_wristfracture
  disutility_wristfracture_beta <- (disutility_wristfracture_alpha/disutility_wristfracture) - disutility_wristfracture_alpha
  disutility_wristfracture <- rbeta(n = n_samples, shape1 = disutility_wristfracture_alpha, shape2 = disutility_wristfracture_beta)
  disutility_vertebralfracture <-   0.839 - 0.78
  disutility_vertebralfractureSE <- (((  0.839 - (1.1*0.78)) - (  0.839 - (0.9*0.78)))/3.92)
  disutility_vertebralfracture_alpha <- (disutility_vertebralfracture ^ 2 * (1 - disutility_vertebralfracture)/disutility_vertebralfractureSE ^ 2) - disutility_vertebralfracture
  disutility_vertebralfracture_beta <- (disutility_vertebralfracture_alpha/disutility_vertebralfracture) - disutility_vertebralfracture_alpha
  disutility_vertebralfracture <- rbeta(n = n_samples, shape1 = disutility_vertebralfracture_alpha, shape2 = disutility_vertebralfracture_beta)
  disutility_osteoporosis <- (probability_hipfracture * disutility_hipfracture) + (probability_wristfracture * disutility_wristfracture) + (probability_vertebralfracture * disutility_vertebralfracture)
  # alternative: Li et al. 2022 

  disutility_NHL <- runif(n = n_samples, min = 0.036, max = 0.136) 
  
  disutility_GIC <- runif(n = n_samples, min =  0.839-0.74, max = 0.839-0.68) # min from Djalalov et al 2014 min = stage 1-3, max=stage 4 

  
  disutility_biopsy_adults <- rtri(n = n_samples, min = 0, max = 0.005, mode = 0.003)
  disutility_biopsy_children <- rtri(n = n_samples, min = 0, max = 0.010, mode = 0.006)
  
  disutility_biopsy_wait <- (utility_GFD_children - utility_undiagnosedCD_children) * 6/52 
 
  
  
  #############################################################################
  ## Costs ####################################################################
  #############################################################################
  
  
  cost_osteoporosis <- rtri(n=n_samples, 
                            min = 0.75*((110700000/507200)*1.0398), 
                            max = 1.25*((110700000/507200)*1.0398), 
                            mode = (110700000/507200)*1.0398) # 110,700,000 euros in 2019, 507200 people with OSTP in 2020 ((https://www.vzinfo.nl/osteoporose/zorguitgaven))

  
  
  
  cost_IDA <- (52*(2*0.02))+(4*6) # 0.02 per tab, 2 times per week, 4 times the afleverkosten which = 6 eur

  cost_gfp <- rtri(n=n_samples,
                 min = 0.75*(950*(110.39/107.51)), # based on NCV (van Overveld)
                 max = 1.25*(950*(110.39/107.51)),
                 mode = 950*(110.39/107.51))

  cost_undiagnosedCD_adults <- (421*1.0353)*1.16 # these are march 2021 euros, was 421 gbp in 2019
  cost_undiagnosedCD_se_adults <- (3.34*1.0353)*1.16 # was 3.34 gbp in 2019
  cost_undiagnosedCD_alpha_adults <- (cost_undiagnosedCD_adults/cost_undiagnosedCD_se_adults)^2
  cost_undiagnosedCD_beta_adults <- (cost_undiagnosedCD_se_adults^2)/cost_undiagnosedCD_adults
  cost_undiagnosedCD_adults <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha_adults, scale = cost_undiagnosedCD_beta_adults) 
  
  cost_undiagnosedCD_children <- (248*1.0353)*1.16 # was 248 gbp in 2019
  cost_undiagnosedCD_se_children <- (4.97*1.0353)*1.16
  cost_undiagnosedCD_alpha_children <- (cost_undiagnosedCD_children/cost_undiagnosedCD_se_children)^2
  cost_undiagnosedCD_beta_children <- (cost_undiagnosedCD_se_children^2)/cost_undiagnosedCD_children
  cost_undiagnosedCD_children <- rgamma(n = n_samples, shape = cost_undiagnosedCD_alpha_children, scale = cost_undiagnosedCD_beta_children) 
  
  cost_CDGFD_adults <- (757*1.0353)*1.16
  cost_CDGFD_se_adults <- (5.3*1.0353)*1.16
  cost_CDGFD_alpha_adults <- (cost_CDGFD_adults/cost_CDGFD_se_adults)^2
  cost_CDGFD_beta_adults <- (cost_CDGFD_se_adults^2)/cost_CDGFD_adults
  cost_CDGFD_adults <- rgamma(n = n_samples, shape = cost_CDGFD_alpha_adults, scale = cost_CDGFD_beta_adults)

  cost_CDGFD_children <- (452*1.0353)*1.16
  cost_CDGFD_se_children <- (20.6*1.0353)*1.16
  cost_CDGFD_alpha_children <- (cost_CDGFD_children/cost_CDGFD_se_children)^2
  cost_CDGFD_beta_children <- (cost_CDGFD_se_children^2)/cost_CDGFD_children
  cost_CDGFD_children <- rgamma(n = n_samples, shape = cost_CDGFD_alpha_children, scale = cost_CDGFD_beta_children)

  
  
  
  probability_biopsy_adults <- runif(n = n_samples, min = 0.6, max = 0.8)
  probability_biopsy_children <- rbeta(n=n_samples, shape1 = 151, shape2 = (480 - 151)) # 480=total in NSCK, 151=n who received biopsy
  
  
  
  cost_NHL <- rtri( 
    n = n_samples,
    min = 0.75*((207000000*1.0398)/(12880 +	14127)), # updated based on 20-year prevalence in 2019 and KvZ total in 2019. Was in Bristol: 18396
    max = 1.25*((207000000*1.0398)/(12880 +	14127)),
    mode= (207000000*1.0398)/(12880 +	14127)
  )
  

  cost_diagnosis_ns <- 651.88 #2021 Euros
  cost_diagnosis_sd_ns <- 199.73
  cost_diagnosis_location_ns <- log(cost_diagnosis_ns^2 / sqrt(cost_diagnosis_sd_ns^2 + cost_diagnosis_ns^2))
  cost_diagnosis_shape_ns <- sqrt(log(1 + (cost_diagnosis_sd_ns^2 / cost_diagnosis_ns^2)))
  cost_diagnosis_ns <- rlnorm(n = n_samples, cost_diagnosis_location_ns,  cost_diagnosis_shape_ns)

  
  
  ferritine <- 		6.48 #all below are 2021 euros
  ijzer <-		2.35
  transferrine <-		4.41
  tpo_antistoffen <-		31.71
  vit_b12 <-		6.45
  foliumzuur <-		5.86
  endo_iga <- 		34.96
  iga_anti_ttg <-		34.96
  hla <- 		263.44
  order_tarif <-  12.13
  biopsy <- 376.65
  
  cost_diagnosis_s <-   ferritine + ijzer + transferrine +tpo_antistoffen +  vit_b12 +foliumzuur +endo_iga +iga_anti_ttg +hla + order_tarif + (biopsy*(5/59)) # 5/59 = prelim. prop receiving biopsy in POCT
  
  
  
  cost_GIC <- rtri(
   n = n_samples,
   min = 0.75*((558000000*1.0398)/108465), #based on 20-year prevalence in 2019 and total KvZ for colon cancer in 2019. Links: https://www.vzinfo.nl/dikkedarmkanker/zorguitgaven & https://iknl.nl/nkr-cijfers?fs%7Cepidemiologie_id=528&fs%7Ctumor_id=138&fs%7Cprevalentie_id=554&fs%7Cperiode_id=563&fs%7Cgeslacht_id=644&fs%7Cleeftijdsgroep_id=677&fs%7Cjaren_na_diagnose_id=687&cs%7Ctype=false&cs%7CxAxis=false&cs%7Cseries=epidemiologie_id&ts%7CrowDimensions=&ts%7CcolumnDimensions=&lang%7Clanguage=nl
   max = 1.25*((558000000*1.0398)/108465),
   mode = ((558000000*1.0398)/108465)
  )
  
  


sens_POCT <- logit(0.98)
sens_POCT_lci <- logit(0.94)
sens_POCT_uci <- logit(0.99)
sd = (sens_POCT_uci-sens_POCT_lci)/3.92
sens_POCT <- rnorm(n=n_samples, sens_POCT, sd=sd)
sens_POCT <- expit(sens_POCT)

spec_POCT <- logit(0.97)
spec_POCT_lci <- logit(0.94)
spec_POCT_uci <- logit(0.99)
sd = (spec_POCT_uci-spec_POCT_lci)/3.92
spec_POCT <- rnorm(n=n_samples, spec_POCT, sd=sd)
spec_POCT <- expit(spec_POCT)

  
  cost_POCT <- 11.51 + 13 # the 13 is for the POCT itself, the rest is for related costs. 2021 prices
  
  
  pre_test_probability_overall <- rtri(
                                        n = n_samples,
                                        min = 0.008,
                                        max = 0.015,
                                        mode = 0.01) 
  
  cost_quest <- 2 # cost of questionnaire per child
  
  #######################
  input_parameters <- data.frame(probability_late_diagnosis_children, probability_late_diagnosis_adults, 
                                 probability_osteoporosis, probability_NHL, probability_GIC, probability_nocomplications, probability_IDA,
                                 osteoporosis_lograte, log_or_osteoporosis_GFD, log_or_osteoporosis_noGFD,
                                 NHL_lograte, log_rr_NHL_GFD, log_rr_NHL_noGFD,
                                 GIC_lograte, log_rr_GIC_GFD, log_rr_GIC_noGFD,
                                 death_probability_NHL, death_log_hr_osteoporosis_mixedj, death_probability_GIC, 
                                 utility_GFD_children, utility_GFD_adults, utility_undiagnosedCD_children, utility_undiagnosedCD_adults,  disutility_osteoporosis, disutility_NHL, disutility_GIC, #disutility_fp,
                                 probability_biopsy_children, probability_biopsy_adults, disutility_biopsy_adults, disutility_biopsy_children, disutility_biopsy_wait,
                                 cost_CDGFD_children, cost_CDGFD_adults, cost_osteoporosis, cost_undiagnosedCD_children, cost_undiagnosedCD_adults, cost_IDA, cost_GIC, cost_NHL,
                                 cost_diagnosis_ns, cost_diagnosis_s, #ns= no screening, s=screening  
                                 sens_POCT, spec_POCT, cost_POCT, cost_gfp, pre_test_probability_overall, cost_quest)
  
  return(input_parameters)
} #END FUNCTION

