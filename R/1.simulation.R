#######################################################################
##                    Simulate household epidemics
##
## author: Maylis Layan
## creation date: 2020/05/01
#######################################################################

rm(list = ls())
setwd("../../")
library(tidyverse)
library(Rcpp)

####################################################
## Prepare environment
####################################################
# Remove any .o file (produced by windows)
filesToRemove = list.files(path = "cpp/simulation/", pattern = "^.*\\.o$", full.names = T)
if (length(filesToRemove)) sapply(filesToRemove, file.remove)

# Source cpp file
sourceCpp("cpp/simulation/simulate_epidemic.cpp")

####################################################
## Load and prepare household data
####################################################
# Load input data
bdd = read.table("dta/2021_05_14_model_data_cpp_full_database_2doses.txt",
                 sep = " ", stringsAsFactors = F, header = F)
colnames(bdd) = c("indid", "hhid", "hhsize", "dds", "infectionStatus", "vaccinated", "studyPeriod", "adult", "index", "isolation")

# Modify index case when dds of contact is anterior to dds of identified index case
bdd$index = sapply(bdd$indid, function(x) {
  h = bdd$hhid[bdd$indid == x]
  sub = bdd[bdd$hhid == h, ]
  dds = sub$dds
  d = sub$dds[sub$indid == x]
  i = out = sub$index[sub$indid == x]
  
  if (sum(sub$dds != 1000 & sub$index == 0) > 0 ) {
    if (d != 1000 & d < min(sub$dds[sub$index == 1]) & i == 0) out = 1
    if (i == 1 & d > min(sub$dds[sub$index == 0 & sub$dds != 1000])) out = 0
  }
  
  return(out)
}) 

# Proportion of asymptomatic cases among household contacts
pAsymptomatic = sapply(
  c(2:8,12), 
  function(x) {
    sum(bdd$infectionStatus == 2 & bdd$index == 0 & bdd$hhsize ==x) /
    sum(bdd$index == 0 & bdd$dds != 1000  & bdd$hhsize ==x)
  }
)

# Convert into list format
hhList = unique(bdd$hhid)
data = lapply(hhList, function(x) bdd[bdd$hhid ==x, ])

# Input data for simulation
bdd$dds[bdd$index == 0] = 1000
bdd$infectionStatus[bdd$index == 0] = 0
inputData = lapply(hhList, function(x) bdd[bdd$hhid == x, ])

# HH with only 1 identified index case
hhids_1_index = sapply(unique(bdd$hhid), function(x) sum(bdd$index[bdd$hhid %in% x]))
hhids_1_index = unique(bdd$hhid)[hhids_1_index == 1]

############################################
## Load posterior distributions
############################################
parameterValue = read.table("results/mcmc_full_database_2doses_1.0_1.0_0.6_1.txt", 
                            header = T)

parameterValue = parameterValue[parameterValue$iteration >= 10000, ]

############################################
## Simulate household clusters using
## parameter values drawn from posteriors
############################################
nSim = 2000
expected_sar = data.frame()
expected_dist = data.frame()

seed <- as.integer(Sys.time() + Sys.getpid())
print(paste('Random seed =', seed))
set.seed(seed)

for (i in 1:nSim) {
  
  # Draw parameters
  k = sample(1:nrow(parameterValue), 1)
  b = parameterValue$beta[k]
  a = parameterValue$alpha[k]
  delta = parameterValue$delta[k]
  rSAV = parameterValue$rSAV[k]
  rSAI = parameterValue$rSAI[k]
  rSAVI = parameterValue$rSAVI[k]
  rSC = parameterValue$rSC[k]
  rSCI = parameterValue$rSCI[k]
  rInfVac = parameterValue$rInfVac[k]
  
  # Simulate epidemic
  out = lapply(inputData, 
               hhEpidemic,
               beta = b,
               alpha = a,
               delta = delta,
               rSAVI = rSAVI,
               rSAV = rSAV,
               rSAI = rSAI,
               rSC = rSC,
               rSCI = rSCI,
               rInfVac = rInfVac,
               pAsymptomatic = pAsymptomatic,
               mainHHSize = 4.0,
               dt = 0.01)
  
  # Store expected secondary attack rate and number of secondary cases
  expected_sar_temp = do.call("rbind", out) %>%
    filter(hhid %in% hhids_1_index) %>%
    group_by(hhid) %>%
    summarise(hhsize = n(), nInf = sum(infectionStatus > 0 & dds <= studyPeriod) - sum(index), nContacts = sum(index == 0)) %>%
    group_by(hhsize) %>%
    summarise(nInfContacts = sum(nInf), nContacts = sum(nContacts), .groups = "drop")
  
  expected_dist_temp = do.call("rbind", out) %>%
    group_by(hhid) %>%
    summarise(hhsize = n(), nInf = sum(infectionStatus > 0 & dds <= studyPeriod)) %>%
    group_by(hhsize, nInf) %>%
    summarise(nHH = n(), .groups = "drop")
  
  expected_sar_temp$sim = expected_dist_temp$sim = i
  expected_sar = bind_rows(expected_sar, expected_sar_temp)
  expected_dist = bind_rows(expected_dist, expected_dist_temp)
  
  if (i %% 100 == 0) print(i)
}

############################################
## Write outputs
############################################
write.csv(expected_sar,
          "results/expectedSAR_full_database_2doses.csv",
          row.names = F)

############################################
## Table comparing observed and expected 
## distributions
############################################
# Observed number of cases per household size
observed_dist = do.call("rbind", data) %>%
  group_by(hhid) %>%
  summarise(hhsize = n(), nInf = sum(infectionStatus > 0)) %>%
  group_by(hhsize, nInf) %>%
  summarise(mean =n(), .groups = "drop") %>%
  filter(hhsize <= 5) %>%
  mutate(distribution = "Observed")

# Comparison
expected_dist %>%
  filter(hhsize <= 5) %>%
  select(-sim) %>%
  group_by(hhsize, nInf) %>%
  summarise(
    mean = round(mean(nHH),1), 
    q025 = round(quantile(nHH, 0.025),1),
    q975 = round(quantile(nHH, 0.975),1)
  ) %>%
  mutate(distribution = "Expected") %>%
  bind_rows(., observed_dist) %>%
  mutate(value = ifelse(is.na(q025), paste(mean), paste0(mean, " (", q025, "-", q975, ")"))) %>%
  select(-c(mean, q025, q975)) %>%
  pivot_wider(names_from = nInf, values_from = value) %>%
  arrange(hhsize, desc(distribution)) %>%
  replace(is.na(.), "") %>%
  write.table(., "results/comparedDistributions_full_database_2doses.csv",
              row.names = F)
