#!/usr/bin/env Rscript

library(dplyr, warn.conflicts = F, quietly = T)

# Set parameters of the mcmc chain
#   Model: 
#     - baseline: estimate the risk of infection in the community and in the households 
#     - rInf: estimate the risk of infection in the community and in the households while accounting for the relative infectivity of vaccinated cases
#     - rSC: estimate the risk of infection in the community and in the households while accounting for the relative susceptibility of the different contact categories 
#     - rInf_rSC: estimate risk of infection in the community and in the households while accounting for the relative infectivity of vaccinated cases and the relative susceptibility of the different contact categories
#   In the paper, we used the rInf_rSC model
#   
#   Frequency-dependent transmission:
#     - 1: frequency-dependent transmission
#     - 0: density-dependent transmission
#   
#   Estimation of the power coefficient of the frequency-dependent transmission formulation:
#     - 1: estimation of delta
#     - 0: set to its default value 1
#   
#   Chain ID:
#     - 1, 2 or 3
#   
#   Maximum PCR detectability:
#     - 10: up to 10 days after infection
#     - 10: up to 15 days after infection
#   
#   Log-sd of the relative infectivity of vaccinated cases compared to the unvaccinated cases:
#     - In the paper, we used the following values: 0.7, 1, 2
#   
#   Log-sd of the relative susceptibility of the different contact categories compared to the unvaccinated and not isolated adult/teenager contacts:
#     - In the paper, we used the following values: 0.7, 1, 2
#   
#   Relative infectivity of asymptomatic cases compared to the symptomatic cases:
#     - In the paper, we used the following values: 0.6, 1
#   
#   Vaccination definition:
#     - 1dose: vaccine is effective from 15 days after the 1st dose
#     - 2doses: vaccine is effective from 7 days after the 2nd dose 
#   
#   Database to be analyzed:
#     - full_database: principal analysis
#     - 1PCR: households where all negative contacts performed at least one PCR test in the 10 days following the detection of the index case
#     - 2PCR: households where all negative contacts performed at least two PCR tests in the 10 days following the detection of the index case
#     - nor_early_vaccination: does not contains the households where the index case was vaccinated but got infected before the vaccine was considered effective (>7 days after the 2nd dose)
#   
#   Main household size:
#     - We scaled the frequency-dependent formulation with the main household size of the database that is analyzed to de-correlate the power coefficient and the baseline risk of infection in the households
  

  
baseline = expand.grid(
  "rS_rInfVac", # Models
  1, # Frequency dependent transmission
  1, # estimation of delta 
  1:3, # Chain ID
  10, # maxPCRDetectability
  1, # sdrInfVac
  1, # sdrSVac
  0.6, # rAsymptomatic
  "2doses", # Vaccination definition
  "full_database", # data base
  4, # main household size
  stringsAsFactors = F
)

# Add sensitivity analysis
params = baseline

for (v in c(0.7,2)) { # Prior log-sd
  for (col in c("Var6", "Var7")) {
    params_temp = baseline
    params_temp[[col]] = v
    params = bind_rows(params, params_temp)
  }
}

# Asymptomatic cases as infectious as symptomatic cases
params_temp = baseline
params_temp$Var8 = 1
params = bind_rows(params, params_temp)

# Vaccine effective from 15 days after 1st dose
params_temp = baseline
params_temp$Var9 = "1dose"
params = bind_rows(params, params_temp)

# 1 PCR tests for all contacts
params_temp = baseline
params_temp$Var11 = "known_outcome"
params_temp$Var12 = 3
params = bind_rows(params, params_temp)

# 2 PCR tests for all negative contacts
params_temp = baseline
params_temp$Var11 = "2PCR"
params_temp$Var12 = 3
params = bind_rows(params, params_temp)

# No index case with early vaccination
params_temp = baseline
params_temp$Var11 = "no_early_vaccination"
params_temp$Var12 = 4
params = bind_rows(params, params_temp)


# Build sh files
for (curr in 1:nrow(params)) {

  # Command to run the script
  arg_string = paste(as.vector(t(params[curr, ])), collapse = " ")
  system(paste("./program.out", arg_string))
  
}
