###################################################
##            SENSIVITITY ANALYSES 
##        
## author: Maylis Layan
## creation date: 2021/05/06
###################################################

rm(list = ls())
setwd("../")

library(tidyverse)
library(gridExtra)

############################################
## Load outputs
############################################
burnIn = 1000

mcmc_posterior = data.frame()
mcmc_LL = data.frame()
mcmc_accept = data.frame()

parameterNames = c("alpha", "beta", "delta", "rSAI", "rSAV", "rSAVI", "rSC", "rSCI", "rInfVac")

# Scenarios
scenarios = data.frame(
  name = c(
    "Baseline", 
    ">15 days after 1st dose",
    "2 PCR tests for all contacts",
    "100% infectivity of asymptomatic cases",
    "relative infectivity prior with log-sd=0.7", 
    "relative infectivity prior with log-sd=2",
    "relative susceptibility prior with log-sd=0.7", 
    "relative susceptibility prior with log-sd=2"
    ),
  fileName = c(
    "results/mcmc_full_database_2doses_1.0_1.0_0.6_1.txt",
    "results/mcmc_full_database_1dose_1.0_1.0_0.6_1.txt",
    "results/mcmc_known_outcome_2doses_1.0_1.0_0.6_1.txt",
    "results/mcmc_full_database_2doses_1.0_1.0_0.6_1.txt",
    "results/mcmc_full_database_2doses_0.7_1.0_0.6_1.txt",
    "results/mcmc_full_database_2doses_2.0_1.0_0.6_1.txt",
    "results/mcmc_full_database_2doses_1.0_0.7_0.6_1.txt",
    "results/mcmc_full_database_2doses_1.0_2.0_0.6_1.txt"
    )
)

# Load outputs
for (i in 1:nrow(scenarios)) {
  
  # File
  f = scenarios$fileName[i]
  out_temp = read.table(f, sep = " ", header = T)
  
  # LogLik
  LL_temp = subset(out_temp, subset = iteration >= burnIn, select = c(iteration, logLik)) 
  LL_temp$Scenario = scenarios$name[i]
  mcmc_LL = bind_rows(mcmc_LL, LL_temp)
  
  # Acceptance rate 
  nAccepted = sapply(out_temp[, colnames(out_temp)[grepl("_", colnames(out_temp))] ], sum)
  params = colnames(out_temp)[grepl("_a", colnames(out_temp))]
  params = gsub("_a", "", params)
  
  accept_temp = sapply(params, function(x) {
    if (nAccepted[paste0(x, "_p")] == 0) {
      return(NA)
    } else {
      out = nAccepted[paste0(x, "_a")] / nAccepted[paste0(x, "_p")]
      names(out) = NULL
      return(out)
    }
  })
  accept_temp = data.frame(rate = accept_temp, parameter = names(accept_temp), row.names = NULL)
  accept_temp$Scenario = scenarios$name[i]
  mcmc_accept = bind_rows(mcmc_accept, accept_temp)
  
  # Acceptance rate 
  notAnalyzed = names(nAccepted[nAccepted == 0])
  notAnalyzed = unique(gsub("_.*", "", notAnalyzed))
  
  posterior_temp = out_temp[out_temp$iteration >= burnIn, c("iteration", params[params != "data"]) ]
  posterior_temp[, notAnalyzed] = NA
  posterior_temp$Scenario = scenarios$name[i]
  mcmc_posterior = bind_rows(mcmc_posterior, posterior_temp) 
}

##################################################
## Parameter estimates
##################################################
# Figures
mcmc_posterior %>%
  pivot_longer(c(rSAV, rSAVI, rSAI, rSC, rSCI), names_to = "susceptibility", values_to = "value") %>%
  group_by(Scenario, susceptibility) %>%
  summarise(median = median(value), 
            q025 = quantile(value, 0.025), 
            q975 = quantile(value, 0.975)
  ) %>% 
  mutate(
    susceptibility = factor(susceptibility, c("rSAI", "rSAV", "rSAVI", "rSC", "rSCI"), c("Isolated\nadult", "Vaccinated\nadult", "Isolated +\nVaccinated\nadult", "Child", "Isolated\nchild")),
    Scenario = factor(Scenario, scenarios$name)
  ) %>%
  filter(!grepl("relative infectivity prior", Scenario)) %>%
  ggplot(., aes(x = susceptibility, y = median, ymin = q025, ymax = q975, col = Scenario, shape = Scenario)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  theme_light() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12)) +
  scale_shape_manual(values = c(16,15,17,4,18,25,19)) +
  scale_color_manual(values = c("black", "gold", "orange", "red", "pink", "cadetblue2")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits= c(0,1.1)) +
  labs(x = "", y = "Relative susceptibility", col = "", shape = "")
ggsave("figures/parameter_estimates_susceptibility_sensitivity.pdf", width = 9.5, height = 4)

mcmc_posterior %>%
  group_by(Scenario) %>%
  summarise(median = median(rInfVac), 
            q025 = quantile(rInfVac, 0.025), 
            q975 = quantile(rInfVac, 0.975)
  ) %>% 
  mutate(
    rInfVac = "Vaccinated cases",
    Scenario = factor(Scenario, scenarios$name)
    ) %>%
  filter(!grepl("relative susceptibility prior", Scenario)) %>%
  ggplot(., aes(x = rInfVac, y = median, ymin = q025, ymax = q975, col = Scenario, shape = Scenario)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  theme_light() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 12)) +
  scale_shape_manual(values = c(16,15,17,4,18,25,19)) +
  scale_color_manual(values = c("black", "gold", "orange", "red", "pink", "cadetblue2")) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.2), limits= c(0,1.1)) +
  labs(x = "", y = "Relative infectivity", col = "", shape = "")
ggsave("figures/parameter_estimates_infectivity_sensitivity.pdf", width = 7, height = 3.75)
  

# Tables 
mcmc_posterior %>%
  pivot_longer(c(rSAV, rSAVI, rSAI, rSC, rSCI), names_to = "Relative susceptibility", values_to = "value") %>%
  group_by(Scenario, `Relative susceptibility`) %>%
  summarise(median = round(median(value), 2), 
            q025 = round(quantile(value, 0.025), 2), 
            q975 = round(quantile(value, 0.975),2)
  ) %>% 
  mutate(
    `Relative susceptibility` = factor(
      `Relative susceptibility`, 
      c("rSAI", "rSAV", "rSAVI", "rSC", "rSCI"), 
      c("Isolated adult", "Vaccinated adult", "Isolated + vaccinated adult", "Child", "Isolated child")),
    Scenario = factor(Scenario, scenarios$name), 
    value = paste0(median, " (", q025, "-", q975, ")")
  ) %>%
  filter(!grepl("relative infectivity prior", Scenario))  %>%
  pivot_wider(-c(median, q025, q975), values_from = value, names_from = Scenario) %>%
  write.table(., "tables/parameter_estimates_susceptibility_sensitivity.csv", sep = ",", row.names = F)
  

mcmc_posterior %>%
  group_by(Scenario) %>%
  summarise(median = round(median(rInfVac),2), 
            q025 = round(quantile(rInfVac, 0.025),2), 
            q975 = round(quantile(rInfVac, 0.975),2)
  ) %>% 
  mutate(
    `Relative infectivity` = "Vaccinated case",
    Scenario = factor(Scenario, scenarios$name),
    value = paste0(median, " (", q025, "-", q975, ")")
  ) %>%
  filter(!grepl("relative susceptibility prior", Scenario)) %>%
  pivot_wider(-c(median, q025, q975), values_from = value, names_from = Scenario) %>%
  write.table(., "tables/parameter_estimates_infectivity_sensitivity.csv", sep = ",", row.names = F)
  
