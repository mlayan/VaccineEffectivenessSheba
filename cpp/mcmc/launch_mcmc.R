# Baseline 
baseline = expand.grid(
  "rS_rInfVac", # In the paper, we use the model that estimates the relative infectivity and relative susceptibility parameters, it is rS_rInfVac
   1, #0: density-dependent transmission; 1: frequency-dependent transmission (used in the manuscript)
   1, #0: default value of delta (power coefficient on household size in the transmission risk) set to 1; 1: delta is estimated (used in the manuscript)
   1, #1; 2; 3: Chain ID
   10, # Maximum period of detectability by PCR for asymptomatic cases, we used 10 days in the manuscript
   1, # log-sd of the relative infectivity of vaccinated cases compared to unvaccinated cases. We used 0.7, 1 (baseline), and 2 in the manuscript.
   1, # log-sd of the relative susceptibility of the different contact categories We used 0.7, 1 (baseline), and 2 in the manuscript.
   0.6, # relative infectivity of asymptomatic cases compared to symptomatic cases, we used 0.6 (baseline) and 1 (sensitivity analysis) in the manuscript
   "2doses", # //1dose: effective vaccination if exposure occurs >=15 days after the 1st dose; 2doses: effective vaccination if exposure occurs >= 7 days after 2nd dose
   "full_database", # full_database; 1PCR; 2PCR; no_early_vaccination: name of the data base 
   4, # main household size in the database that is analyzed
   stringsAsFactors = F
)
params = baseline

# Sensitivity analysis
# Modify log-sd of relative parameters
# params = baseline
# params[[c("Var6", "Var7")]] = 0.7
# params[[c("Var6", "Var7")]] = 2

# Asymptomatic cases are as infectious as symptomatic cases
# params = baseline
# params$Var8 = 1

# Vaccine effective from 15 days after 1st dose
# params= baseline
# params$Var9 = "1dose"

# 1 PCR tests for all contacts
# params_temp = baseline
# params_temp$Var10 = "1PCR"
# params_temp$Var11 = 3

# 2 PCR tests for all negative contacts
# params = baseline
# params_temp$Var10 = "2PCR"
# params_temp$Var11 = 3

# No index case with early vaccination
# params = baseline
# params$Var10 = "no_early_vaccination"
# params$Var11 = 4

# Run the model you want
system(paste0("program.out ", paste0(t(params), collapse = " ")))

