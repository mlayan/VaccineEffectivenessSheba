###################################################
##                  ANALYSIS
##                MODEL ADEQUACY
##
## author: Maylis Layan
###################################################

rm(list = ls())
setwd("../") # Give path of the project

library(tidyverse)
library(gt)
library(gridExtra)
library(binom)


############################################
## Supplementary figure 1 with model fit
############################################
## Load data
data = read.table(paste0("data/2021_05_14_full_database_2doses.txt"), sep = " ", header = F)
colnames(data) = c("indid", "hhid", "hhsize", "dds", "infectionStatus", "vaccinated", "studyPeriod", "adult", "index", "isolation")

# Observed SAR
obsData = data %>%
  filter(index == 0) %>% 
  group_by(hhsize) %>%
  summarise(
    median = sum(infectionStatus > 0) / n(),
    q025 = binom.confint(sum(infectionStatus > 0), n(), methods = "exact")[["lower"]],
    q975 = binom.confint(sum(infectionStatus > 0), n(), methods = "exact")[["upper"]]
    ) %>%
  mutate(data = "Observed") %>%
  filter(hhsize <= 5) %>%
  data.frame()

# Expected SAR
sar_exp = read.table("results/2doses/sar_rS_rInfVac_full_database.txt", header=T) %>%
  mutate(data = "Expected")

ggplot(bind_rows(obsData, sar_exp)) +
  geom_pointrange(aes(x = as.factor(hhsize), y = median, ymin = q025, ymax = q975, col = data), 
                  position = position_dodge(width = 0.5), fatten = 3) +
  theme_light() +
  scale_color_manual(values = c("black", "cornflower blue")) +
  ylim(c(0,1)) +
  labs(x = "Household size", y = "Secondary attack rate", col = "", shape = "")
ggsave(paste0("figures/2doses/model_adequacy.png"), height = 3.5, width = 5)


############################################
## Compare observed and expected distributions
## Supplementary Table 1
############################################
# Observed distribution
observedDist = read.table("data/2021_05_14_full_database_2doses.txt", sep = " ", header = F) %>%
  rename(indid = V1, hhid = V2, hhsize = V3, dds = V4, infectionStatus = V5, vaccinated = V6, studyPeriod = V7, 
         adult = V8, index = V9, isolation = V10) %>%
  group_by(hhid) %>%
  summarise(hhsize = n(), nInf = sum(infectionStatus > 0), .groups = "drop") %>%
  group_by(hhsize, nInf) %>%
  summarise(mean = n(), .groups = "drop")  %>%
  mutate(distribution = "Observed")

# Expected distribution
temp_SAR = read.table("results/2doses/infection_dist_rS_rInfVac_full_database.txt", header = T) %>%
  bind_rows(., observedDist) %>%
  mutate(value = ifelse(is.na(q025), paste(mean), paste0(round(mean, 1), " (", q025, "-", q975, ")"))) %>%
  select(-c(mean, q025, q975)) %>%
  filter(hhsize <= 5) %>%
  pivot_wider(names_from = nInf, values_from = value) %>%
  arrange(hhsize, desc(distribution)) %>%
  replace(is.na(.), "")
    
write.table(temp_SAR, "tables/comparedDistributions_2doses_full_database_rS_rInfVac.csv", 
            row.names = F)
  
  
  