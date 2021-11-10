#######################################################################
##                  MCMC - Transmission model
##
## author: Maylis Layan
#######################################################################

rm(list = ls())

setwd("../") # Give path of the project
library(dplyr)
library(tidyr)
library(Rcpp)

# Source cpp file
sourceCpp("cpp/simulation/simulate_epidemic.cpp")

####################################################
## Load and prepare household data
####################################################
# Load data
bdd = read.table(paste0("data/2021_05_14_full_database_2doses.txt"), 
                 sep = " ", stringsAsFactors = F, header = F)
colnames(bdd) = c("indid", "hhid", "hhsize", "dds", "infectionStatus", "vaccinated", 
                  "studyPeriod", "adult", "index", "isolation")
mainHHSize = 4

# Proportion of asymptomatic cases among household contacts
pAsymptomatic = sum(bdd$infectionStatus == 2 & bdd$index == 0) / 
  sum(bdd$index == 0 & bdd$dds != 1000)
hhComp = table(bdd$hhsize[bdd$index == 1])
hhids = unique(bdd$hhid)

# Rename households 
hhList = unique(bdd$hhid) 

# Convert into list format
data = lapply(hhList, function(x) bdd[bdd$hhid ==x, ])

# Input data for simulation
bdd$dds[bdd$index == 0] = 1000
bdd$infectionStatus[bdd$index == 0] = 0

inputData = lapply(hhList, function(x) {
  
  # Select household 
  out = bdd[bdd$hhid == x, ]
  
  # Order data frame to have index case on the 1st row
  out = out[order(out$index, decreasing = T), ] 
  
  # Study period is the one of the index case
  out$studyPeriod = max(out$studyPeriod)
  
  return(out)
}
)

############################################
## Load posterior distributions
############################################
fileName = paste0("results/2doses/mcmc_full_database_rS_rInfVac_1_10_1.0_1.0_0.6_1.txt")
parameterValue = read.table(fileName, sep = " ", header = T)

# Remove burn-in of the posterior distribution
parameterValue = parameterValue[ceiling(0.1*nrow(parameterValue)):nrow(parameterValue), ]
row.names(parameterValue) = NULL

############################################
## Simulate household clusters using 
## parameter values drawn from posteriors
############################################
nSim = 2000
simulated_dist = data.frame()

seed <- as.integer(Sys.time() + Sys.getpid())
print(paste('Random seed =', seed))
set.seed(seed)

for (i in 1:nSim) {
  
  # Draw parameters in the posterior distribution
  params = parameterValue[sample(1:nrow(parameterValue), 1), ] 
  beta = params[["beta"]]
  alpha = params[["alpha"]]
  delta = params[["delta"]]
  rSAV = params[["rSAV"]]
  rSAVI = params[["rSAVI"]]
  rSAI = params[["rSAI"]]
  rSC = params[["rSC"]]
  rSCI = params[["rSCI"]]
  rInfVac = params[["rInfVac"]]
  rAsymptomatic = params[["rAsymptomatic"]]
  
  # Simulate epidemic
  out = lapply(inputData, 
               hhEpidemic,
               beta = beta,
               alpha = alpha,
               delta = delta, 
               rSAV = rSAV,
               rSAVI = rSAVI,
               rSAI = rSAI,
               rSC = rSC,
               rSCI = rSCI,
               rInfVac = rInfVac,
               rAsymptomatic = rAsymptomatic,
               pAsymptomatic = pAsymptomatic,
               mainHHSize = mainHHSize, 
               dt = 0.01)
  
  # Store expected number of cases
  sim_i = do.call("rbind", out)
  sim_i$infectionStatus[sim_i$dds > sim_i$studyPeriod & sim_i$dds != 1000] = 0
  sim_i$sim = i
  simulated_dist = bind_rows(simulated_dist, sim_i)
  
  if (i %% 100 == 0) print(i)
}


# SAR per household size
simulated_dist %>%
  filter(hhsize <= 5) %>%
  group_by(sim, hhsize) %>%
  summarise(sar = sum(infectionStatus > 0 & index == 0) / sum(index == 0)) %>%
  group_by(hhsize) %>%
  summarise(
    min = min(sar),
    mean = mean(sar), 
    median = median(sar), 
    q025 = quantile(sar, 0.025), 
    q975 = quantile(sar, 0.975), 
    max = max(sar)
  ) %>%
  write.table(paste0("results/2doses/sar_rS_rInfVac_full_database.txt"), 
              sep = " ", row.names = F, fileEncoding = "UTF-8")


# Distribution of the number of infected individuals per 
# household size
simulated_dist %>%
  filter(hhsize <= 5) %>%
  group_by(sim, hhid) %>%
  summarise(hhsize = n(), nInf = sum(infectionStatus > 0)) %>%
  group_by(sim, hhsize, nInf) %>%
  summarise(nind = n()) %>%
  group_by(hhsize, nInf) %>%
  summarise(
    mean = mean(nind),
    q025 = quantile(nind, 0.025), 
    q975 = quantile(nind, 0.975)
    ) %>%
  mutate(distribution = "Expected") %>%
  write.table(paste0("results/2doses/infection_dist_rS_rInfVac_full_database.txt"), 
              sep = " ", row.names = F, fileEncoding = "UTF-8")
  



