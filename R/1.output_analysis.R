###################################################
##                  ANALYSIS
##            BASIC STATISTICS
##
## author: Maylis Layan
###################################################

rm(list = ls())
setwd("../") # Give path of the project

library(tidyverse)
library(gt)
library(gridExtra)


############################################
## Load MCMC outputs
############################################
# Remove 10% of burn in 
burnIn = 0.1

mcmc_posterior = data.frame()
mcmc_LL = data.frame()
mcmc_accept = data.frame()

parameterNames = c("alpha", "beta", "delta", "rSAI", "rSAV", "rSAVI", "rSC", "rSCI", "rInfVac")

for (vacDef in c("1dose", "2doses")) {
  for (f in list.files(paste0("results/", vacDef), pattern = "^mcmc_.*_10_1.0_1.0_0.6_[0-9].txt$", full.names = T)) {
    
    # Get run info
    if (grepl("PCR", f)) {
      if (grepl('2PCR', f)) database = "2PCR"
      if (grepl('1PCR', f)) database = "1PCR"
      model = gsub("Results.*/mcmc_[0-9]PCR_|_[0-9]{1}.*", "", f)
      asympPeriod = gsub("Results.*/mcmc_[0-9]PCR_[a-zA-Z_0-9]+_|_[0-9]{1}\\.[0-9]{1}_.*.txt", "", f)
      frequencyDependent = gsub("Results.*/mcmc_[0-9]PCR_[a-zA-Z_]+_|_[0-9]{2}.*.txt", "", f)
      sdLnorm = as.numeric(gsub("Results.*/mcmc_[0-9]PCR_[a-zA-Z_]+_[0-9]{1}_[0-9]{2}_|_[0-9]\\.[0-9]_[0-9]\\.[0-9]_[0-9].txt", "", f))
    } else {
      database = regmatches(f, gregexpr("(?<=mcmc_)[a-z_]+(?=_)", f, perl = TRUE))[[1]]
      model = gsub("Results.*/mcmc_[a-z_]+_|_[0-9]{1}.*", "", f)
      asympPeriod = gsub("Results.*/mcmc_[a-zA-Z_]+_[0-9]+_|_[0-9]{1}\\.[0-9]{1}_.*.txt", "", f)
      frequencyDependent = gsub("Results.*/mcmc_[a-zA-Z_]+_|_[0-9]{2}.*.txt", "", f)
      sdLnorm = as.numeric(gsub("Results.*/mcmc_[a-zA-Z_]+_[0-9]{1}_[0-9]{2}_|_[0-9]\\.[0-9]_[0-9]\\.[0-9]_[0-9].txt", "", f))
    }
    chain = gsub("Results.*/.*[0-9]\\.[0-9]_|.txt", "", f)
    
    # Load output data
    out_temp = read.table(f, sep = " ", header = T)
    
    # LogLik
    b=ceiling(burnIn * nrow(out_temp))
    LL_temp = out_temp %>%
      select(iteration, logLik) %>%
      slice(b:n())
    LL_temp$chain = chain
    LL_temp$model = model
    LL_temp$sdLnorm = sdLnorm
    LL_temp$hhSize = frequencyDependent
    LL_temp$vacDef = vacDef
    LL_temp$database = database
    LL_temp$asympPeriod = asympPeriod
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
    accept_temp$chain = chain
    accept_temp$sdLnorm = sdLnorm
    accept_temp$model = model
    accept_temp$hhSize = frequencyDependent
    accept_temp$vacDef = vacDef
    accept_temp$asympPeriod = asympPeriod
    accept_temp$database = database
    mcmc_accept = bind_rows(mcmc_accept, accept_temp)
    
    # Acceptance rate 
    notAnalyzed = names(nAccepted[nAccepted == 0])
    notAnalyzed = unique(gsub("_.*", "", notAnalyzed))
    
    posterior_temp = out_temp %>%
      slice(b:n()) %>%
      select(iteration, all_of(parameterNames))
    posterior_temp[, notAnalyzed] = NA
    posterior_temp$chain = chain
    posterior_temp$sdLnorm = sdLnorm
    posterior_temp$model = model
    posterior_temp$hhSize = frequencyDependent
    posterior_temp$vacDef = vacDef
    posterior_temp$database = database
    posterior_temp$asympPeriod = asympPeriod
    mcmc_posterior = bind_rows(mcmc_posterior, posterior_temp)
    
  }  
}

############################################
## Acceptance rates for each parameter
############################################
mcmc_accept = mcmc_accept %>%
  mutate(parameter = factor(parameter, c(parameterNames, "data")))

ggplot(mcmc_accept[!is.na(mcmc_accept$rate), ],
       aes(x = parameter, y = rate, col = chain)) +
  facet_grid(database ~ vacDef) +
  geom_point(position = position_dodge(width = 0.5)) +
  ylim(c(0,1)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "", y = "Acceptance rate")
ggsave("figures/acceptance_rate_main_analysis.pdf", height = 5, width = 8)


############################################
## Visual inspection of the mixing
############################################
params = params[1:9]

for (v in c("2doses", "1dose")){
  for (d in c("full_database", "1PCR", "2PCR", "no_early_vaccination")) {
    
    k = mcmc_LL$iteration[mcmc_LL$vacDef == v & mcmc_LL$database == d]
    
    if (length(k)) {
      pdf(paste0("figures/", v, "/mixing_", d, "_baseline.pdf"), width = 18, height = 8)
      par(mfrow = c(3, length(params) + 1), mar = c(4,2,3,2))
      for (r in 1:3) {
        plot(
          mcmc_LL$iteration[mcmc_LL$chain == r & mcmc_LL$vacDef == v & mcmc_LL$database == d],
          mcmc_LL$logLik[mcmc_LL$chain == r & mcmc_LL$vacDef == v & mcmc_LL$database == d],
          type = "l",
          xlab = "iteration",
          ylab = "",
          main = "LogLik"
        )
        
        for (p in params) {
          plot(
            mcmc_posterior$iteration[mcmc_posterior$chain == r & mcmc_posterior$vacDef == v & mcmc_posterior$database == d],
            mcmc_posterior[[p]][mcmc_posterior$chain == r & mcmc_posterior$vacDef == v & mcmc_posterior$database == d],
            type = "l",
            xlab = "iteration",
            ylab = "",
            main = p
          )
        }
      }
      dev.off()
    }
  
  }
}


############################################
## Impact of isolation in vaccinated 
## adult/teenager contacts
############################################
# Baseline 
mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "full_database", vacDef == "2doses") %>%
  summarise(
    p = sum(rSAVI < rSAV)/n(),                    # Probability that vaccinated and isolated adult/teenager contacts are less susceptible than vaccinated and not isolated adult/teenager contacts
    median = round(median(rSAVI/rSAV), 2),        # Median relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    q025 = round(quantile(rSAVI/rSAV, 0.025), 2), # Lower bound of the relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    q975 = round(quantile(rSAVI/rSAV, 0.975), 2)  # Upper bound of the relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    ) %>%
  mutate(
    bf = p/(1-p),                                 # Bayes factor 
    evidence = log10(p/(1-p))                     # Level of evidence
  )

# At least 1 PCR for all negative contacts
mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "1PCR", vacDef == "2doses") %>%
  summarise(
    p = sum(rSAVI < rSAV)/n(),                    # Probability that vaccinated and isolated adult/teenager contacts are less susceptible than vaccinated and not isolated adult/teenager contacts
    median = round(median(rSAVI/rSAV), 2),        # Median relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    q025 = round(quantile(rSAVI/rSAV, 0.025), 2), # Lower bound of the relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    q975 = round(quantile(rSAVI/rSAV, 0.975), 2)  # Upper bound of the relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
  ) %>%
  mutate(
    bf = p/(1-p),                                 # Bayes factor 
    evidence = log10(p/(1-p))                     # Level of evidence
  )

# At least 2 PCRs for all negative contacts
mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "2PCR", vacDef == "2doses") %>%
  summarise(
    p = sum(rSAVI < rSAV)/n(),                    # Probability that vaccinated and isolated adult/teenager contacts are less susceptible than vaccinated and not isolated adult/teenager contacts
    median = round(median(rSAVI/rSAV), 2),        # Median relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    q025 = round(quantile(rSAVI/rSAV, 0.025), 2), # Lower bound of the relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
    q975 = round(quantile(rSAVI/rSAV, 0.975), 2)  # Upper bound of the relative susceptibility of vaccinated and isolated adult/teenager contacts compared to vaccinated and not isolated adult/teenager contacts
  ) %>%
  mutate(
    bf = p/(1-p),                                 # Bayes factor 
    evidence = log10(p/(1-p))                     # Level of evidence
  )


############################################
## Fig 2 - 2 doses - Full database
############################################
# Relative susceptibility of the different contact categories
baseline = data.frame(
  rS = "Adult",
  mean = 1,
  median = 1, 
  q025 = 1, 
  q975 = 1
)

G1 = mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", vacDef == "2doses", database == "full_database") %>%
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
                c("Isolated-\nVaccinated-\nAdult", "Isolated+\nVaccinated-\nAdult", 
                  "Isolated-\nVaccinated+\nAdult", "Isolated+\nVaccinated+\nAdult", 
                  "Isolated-\n-\nChild", "Isolated+\n-\nChild"))
  ) %>% 
  ggplot(., aes(x = rS, y = median, ymin = q025, ymax = q975)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_pointrange(position = position_dodge(0.6), col = "#358597", size = 0.8, fatten = 2) + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 9)) + 
  scale_color_manual(values = c("#358597")) +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1.1,0.2)) + 
  labs(y = "", x = "", col = "", title = "Relative susceptibility of contacts", linetype = "")


# Relative infectivity of the different case categories
baseline = data.frame(
  type = "Unvaccinated\n\n",
  mean = 1, 
  median = 1, 
  q025 = 1, 
  q975 = 1
)

G2 = mcmc_posterior %>% 
  filter(chain == "1", model == "rS_rInfVac", database == "full_database", vacDef == "2doses") %>%
  summarise(mean = mean(rInfVac), 
            median = median(rInfVac), 
            q025 = quantile(rInfVac, 0.025), 
            q975 = quantile(rInfVac, 0.975)
  ) %>%
  mutate(type = "Vaccinated\n\n") %>%
  bind_rows(., baseline) %>%
  ggplot(., aes(x = type, y = median, ymin = q025, ymax = q975)) +
  geom_hline(yintercept = 1, col = "black") +
  geom_pointrange(position = position_dodge(0.6), col= "#358597", size = 0.8, fatten = 2) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), 
        axis.text = element_text(size = 9)) + 
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1,0.2)) +
  labs(y = "", x = "", title = "Relative infectivity of cases")

g = grid.arrange(G1,G2, widths = c(2, 1))
ggsave("figures/2doses/fig2.pdf", g, width = 7.5, height = 4)



# Tables 
tab_sus = mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", vacDef == "2doses", database == "full_database") %>%
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
                     c("Isolated- Vaccinated- Adult", "Isolated+ Vaccinated- Adult", 
                       "Isolated- Vaccinated+ Adult", "Isolated+ Vaccinated+ Adult", 
                       "Isolated- Child", "Isolated+ Child"))
    ) %>%
  select(Category, Parameter)


tab_inf = mcmc_posterior %>% 
  filter(chain == "1", model == "rS_rInfVac", database == "full_database", vacDef == "2doses") %>%
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
  write.table(., "tables/parameter_estimates_baseline.csv", sep=",", col.names = F, row.names = F)

############################################
## Fig 3 - 2 doses - full database
############################################
pInf_input = mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "full_database", vacDef == "2doses") %>%
  select(beta, delta, rSAI, rSAVI, rSAV, rSC, rSCI, rInfVac)

# Compute person-to-person transmission probability for
# vaccinated and unvaccinated cases
# Helper function
p_to_p_transmission <- function(df, case, contact) {
  
  df$B <- with(df, beta)
  
  if (contact != "rSA") df$B <- with(df, B*eval(parse(text=contact)))
  if (grepl("Vaccinated adult", case)) df$B <- with(df, B*rInfVac)
  
  df
}

# Get probabilities for households of size 4 
pInf = data.frame()
for (case in c("Vaccinated adult case", "Unvaccinated case")) {
  for (contact in c("rSA", "rSAI", "rSAVI", "rSAV", "rSC", "rSCI")) {
    
    pInf_temp = p_to_p_transmission(pInf_input, case, contact) %>%
      summarise(median = median(1-exp(-B)),
                q025 = quantile(1-exp(-B), 0.025),
                q975 = quantile(1-exp(-B), 0.975)) %>%
      mutate(case = case, 
             contact = contact)
    pInf = bind_rows(pInf, pInf_temp)
  }
}

# Table with estimates
pInf %>%
  mutate(
    median = round(median * 100, 1), 
    q025 = round(q025 *100, 1), 
    q975 = round(q975*100, 1), 
    contact = factor(contact, 
                     levels = c("rSA", "rSC", "rSCI", "rSAI", "rSAV", "rSAVI"), 
                     labels = c("Isolated- Vaccinated- Adult", "Isolated- Child", "Isolated+ Child", 
                                "Isolated+ Vaccinated- Adult", "Isolated- Vaccinated+ Adult", 
                                "Isolated+ Vaccinated+ Adult"))) %>%
  mutate(value = paste0(median, " (", q025, "-", q975, ")")) %>%
  pivot_wider(-c(median, q025, q975), values_from = value, names_from = case) %>%
  write.table(., "tables/transmission_probability_baseline.csv", sep=",", row.names = F)

# Figure with estimates
pInf %>%
  mutate(contact = factor(contact, 
                          levels = c("rSA", "rSC", "rSCI", "rSAI", "rSAV", "rSAVI"), 
                          labels = c("Isolated-\nVaccinated-\nAdult", "Isolated-\n-\nChild", "Isolated+\n-\nChild", 
                                     "Isolated+\nVaccinated-\nAdult", "Isolated-\nVaccinated+\nAdult", 
                                     "Isolated+\nVaccinated+\nAdult"))) %>%
  ggplot(., aes(x = contact, y = median, ymin = q025, ymax = q975, col = case)) +
  geom_pointrange(position = position_dodge(width = 0.4), size = 0.8, fatten = 2) + 
  theme_light() +
  theme(axis.text.x = element_text(size = 9), 
        axis.title = element_text(size =14),
        axis.text.y = element_text(size = 9), 
        plot.title = element_text(size = 12, hjust = 0.5), 
        legend.text = element_text(size = 9)) + 
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#5B0E2D", "#FFA781")) +
  labs(title = "Probability of person-to-person transmission", y = "", x = "", col = "") 
ggsave("figures/2doses/fig3_p_to_p.pdf", width = 6.5, height = 4.5)
