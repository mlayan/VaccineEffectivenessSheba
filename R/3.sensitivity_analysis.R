###################################################
##            SENSIVITITY ANALYSES 
##        
##
## author: Maylis Layan
## creation date: 2021/05/06
## last update: 
###################################################

rm(list = ls())
# setwd("/mnt/gaia/MMMICovid/maylayan/Israel/")
setwd("V:/maylayan/Israel/")

library(tidyverse)
library(gt)
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
    "Results/2doses/test3/mcmc_full_database_rS_rInfVac_1_10_1.0_1.0_0.6_1.txt",
    "Results/1dose/test3/mcmc_full_database_rS_rInfVac_1_10_1.0_1.0_0.6_1.txt",
    "Results/2doses/test3/mcmc_known_outcome_rS_rInfVac_1_10_1.0_1.0_0.6_1.txt",
    "Results/2doses/test3/mcmc_full_database_rS_rInfVac_1_10_1.0_1.0_1.0_1.txt",
    "Results/2doses/test3/mcmc_full_database_rS_rInfVac_1_10_0.7_1.0_0.6_1.txt",
    "Results/2doses/test3/mcmc_full_database_rS_rInfVac_1_10_2.0_1.0_0.6_1.txt",
    "Results/2doses/test3/mcmc_full_database_rS_rInfVac_1_10_1.0_0.7_0.6_1.txt",
    "Results/2doses/test3/mcmc_full_database_rS_rInfVac_1_10_1.0_2.0_0.6_1.txt"
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
## 95% density intervals of prior distributions
##################################################
sds = c(seq(0.5,1, 0.1) , 2)
  
data.frame(
  sd = sds, 
  q025 = round(qlnorm(0.025, 0, sds), 2), 
  q975 = round(qlnorm(0.975,0,sds), 2)
  ) %>% 
  write.table(., "Paper/sd_95IC.csv", sep = ",", row.names = F)

##################################################
## First version
##################################################
# Susceptibility
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
ggsave("Paper/Fig4/fig4_susceptibility.pdf", width = 9.5, height = 4)


# Infectivity
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
ggsave("Paper/Fig4/fig4_infectivity.pdf", width = 7, height = 3.75)
  

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
  write.table(., "Paper/Supplementary_Tables/susceptibility_estimates_sensitivity_analysis.csv", sep = ",", row.names = F)
  

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
  write.table(., "Paper/Supplementary_Tables/susceptibility_estimates_infectivity_analysis.csv", sep = ",", row.names = F)


##################################################
## Early vaccination
##################################################
# Load output data
posterior = read.table("Results/2doses/strict/mcmc_full_database_rS_rInfVac_1_10_1.0_1.0_0.6_1.txt",
                       sep = " ", header = T) %>%
  slice(1000:n()) 

# Susceptibility
baseline = data.frame(
  rS = "Adult",
  mean = 1,
  median = 1, 
  q025 = 1, 
  q975 = 1
)

G1 = posterior %>%
  select(rSAI, rSAVI, rSAV, rSC, rSCI) %>%
  pivot_longer(c(rSAI, rSAVI, rSAV, rSC, rSCI), names_to = "rS", values_to = "value") %>%
  group_by(rS) %>%
  summarise(mean = mean(value),
            median = median(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975)) %>%
  bind_rows(., baseline) %>%
  mutate(
    rS = factor(rS, 
                c("Adult", "rSAI", "rSAV", "rSAVI", "rSC", "rSCI"), 
                c("Adult", "Isolated\nadult", "Vaccinated\nadult", "Isolated +\nVaccinated\nadult", "Child", "Isolated\nchild"))
  ) %>% 
  ggplot(., aes(x = rS, y = median, ymin = q025, ymax = q975)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_pointrange(position = position_dodge(0.6), col = "#358597") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), 
        axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("#358597")) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1.1)) +
  labs(y = "", x = "", col = "", title = "Relative susceptibility", linetype = "")

# Relative infectivity
baseline = data.frame(
  type = "Unvaccinated\n\n",
  mean = 1, 
  median = 1, 
  q025 = 1, 
  q975 = 1
)

G2 = mcmc_posterior %>% 
  summarise(mean = mean(rInfVac), 
            median = median(rInfVac), 
            q025 = quantile(rInfVac, 0.025), 
            q975 = quantile(rInfVac, 0.975)
  ) %>%
  mutate(type = "Vaccinated\n\n") %>%
  bind_rows(., baseline) %>%
  ggplot(., aes(x = type, y = median, ymin = q025, ymax = q975)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_pointrange(position = position_dodge(0.6), col= "#358597") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1.1)) +
  labs(y = "", x = "", title = "Relative infectivity")

g = grid.arrange(G1,G2, widths = c(2, 1))
ggsave("Paper/Supplementary_Figures/fig2_no_early_vaccinated.pdf", g, width = 7.5, height = 4)

# Tables 
tab_sus = posterior %>%
  select(rSAI, rSAVI, rSAV, rSC, rSCI) %>%
  pivot_longer(c(rSAI, rSAVI, rSAV, rSC, rSCI), names_to = "Category", values_to = "value") %>%
  group_by(Category) %>%
  summarise(median = round(median(value), 2),
            q025 = round(quantile(value, 0.025), 2),
            q975 = round(quantile(value, 0.975),2)) %>%
  mutate(
    Parameter = paste0(median, " (", q025, "-", q975, ")"), 
    Category = factor(Category, 
                      c("Adult", "rSAI", "rSAV", "rSAVI", "rSC", "rSCI"), 
                      c("Adult", "Isolated adult", "Vaccinated adult", "Isolated + vaccinated adult", "Child", "Isolated child"))
  ) %>%
  select(Category, Parameter)


tab_inf = mcmc_posterior %>% 
  summarise(median = round(median(rInfVac), 2), 
            q025 = round(quantile(rInfVac, 0.025), 2), 
            q975 = round(quantile(rInfVac, 0.975),2)
  ) %>%
  mutate(Category = "Vaccinated", 
         Parameter = paste0(median, ' (', q025, "-", q975, ")")) %>%
  select(-c(median, q025, q975))

rbind(
  data.frame(Category = "Contact", Parameter = "Relative susceptibility - median (95% CI)"), 
  tab_sus, 
  data.frame(Category = "Case", Parameter = "Relative infectivity - median (95% CI)"), 
  tab_inf
) %>%
  write.table(., "Paper/Supplementary_Tables/parameter_estimates_no_early_vaccination.csv", sep=",", col.names = F, row.names = F)
  
