###################################################
##        ANALYSIS - BASELINE SCENARIO
##
## author: Maylis Layan
## creation date: 2021/04/14
###################################################

rm(list = ls())
setwd("../")

require(tidyverse)
require(gridExtra)

####################################################
## Load original data
####################################################
bdd = read.table("data/2021_05_14_model_data_cpp_full_database_2doses.txt",
                 sep = " ", stringsAsFactors = F, header = F)
colnames(bdd) = c("indid", "hhid", "hhsize", "dds", "infectionStatus", "vaccinated", "studyPeriod", "adult", "index", "isolation")

############################################
## Load outputs
############################################
burnIn = 1000

mcmc_posterior = data.frame()
mcmc_LL = data.frame()
mcmc_accept = data.frame()

parameterNames = c("alpha", "beta", "delta", "rSAI", "rSAV", "rSAVI", "rSC", "rSCI", "rInfVac")


for (f in list.files("results/", pattern = "^mcmc_full_database_2doses_1.0_1.0_0.6_[0-9].txt$", full.names = T)) {

  # Get run info
  chain = gsub("results.*/.*_|.txt", "", f)

  # Load output data
  out_temp = read.table(f, sep = " ", header = T)

  # LogLik
  LL_temp = subset(out_temp, subset = iteration >= burnIn, select = c(iteration, logLik)) 
  LL_temp$chain = chain
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
  mcmc_accept = bind_rows(mcmc_accept, accept_temp)

  # Acceptance rate 
  notAnalyzed = names(nAccepted[nAccepted == 0])
  notAnalyzed = unique(gsub("_.*", "", notAnalyzed))

  posterior_temp = out_temp[out_temp$iteration >= burnIn, c("iteration", params[params != "data"]) ]
  posterior_temp[, notAnalyzed] = NA
  posterior_temp$chain = chain
  mcmc_posterior = bind_rows(mcmc_posterior, posterior_temp)

}  


############################################
## Acceptance rates
############################################
mcmc_accept = mcmc_accept %>%
  filter(!is.na(rate)) %>%
  mutate(parameter = factor(parameter, c(parameterNames, "data")))

ggplot(mcmc_accept, aes(x = parameter, y = rate, col = chain)) +
geom_jitter() +
  ylim(c(0,1)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "", y = "Acceptance rate")
ggsave("figures/acceptance_rate_baseline.pdf", height = 4, width = 5)

############################################
## Mixing
############################################
pdf("figures/mixing_baseline.pdf", width = 20, height = 8)
par(mfrow = c(3, length(parameterNames) + 1), mar = c(4,2,3,2))
for (r in 1:3) {
  plot(
    mcmc_LL$iteration[mcmc_LL$chain == r],
    mcmc_LL$logLik[mcmc_LL$chain == r],
    type = "l",
    xlab = "iteration",
    ylab = "",
    main = "LogLik"
  )

  for (p in parameterNames) {
    plot(
      mcmc_posterior$iteration[mcmc_posterior$chain == r],
      mcmc_posterior[[p]][mcmc_posterior$chain == r],
      type = "l",
      xlab = "iteration",
      ylab = "",
      main = p
    )
  }
}
dev.off()

############################################
## Verify concordance between chains
############################################
logNormalPrior = data.frame(d = dlnorm(seq(0, 10, 0.01), 1), x=seq(0, 10, 0.01))

ggplot(mcmc_posterior, aes(x = alpha, color = chain)) +
  geom_density() +
  labs(x = "Instantaneous risk of infection in the community", y = "Density") +
  theme_light()
ggsave("figures/alpha_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior, aes(x = beta, color = chain)) +
  geom_density() +
  labs(x = "Person-to-person transmission risk", y = "Density") +
  theme_light()
ggsave("figures/beta_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$delta), ], aes(x = delta, color = chain)) +
  geom_density() +
  labs(x = "Power coefficient of the frequency-dependent parametrization", y = "Density") +
  theme_light()
ggsave("figures/delta_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$rInfVac), ], aes(x = rInfVac, color = chain)) +
  geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rInfVac, na.rm = T), ],
            aes(x = x, y = d, color = "grey")) +
  geom_density() +
  theme_light() + 
  labs(x = "Relative infectivity of vaccinated cases", y = "Density") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey"), labels = c("1", "2", "3", "prior"))
ggsave("figures/rInfVac_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSAV), ], aes(x = rSAV, color = chain)) +
  geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSAV, na.rm = T), ],
            aes(x = x, y = d, color = "prior")) +
  geom_density() +
  theme_light() +
  labs(x = "Relative susceptibility of vaccinated and not isolated adults", y = "Density") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey"), labels = c("1", "2", "3", "prior"))
ggsave("figures/rSAV_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSAVI), ], aes(x = rSAVI, color = chain)) +
  geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSAVI, na.rm = T), ],
            aes(x = x, y = d, color = "prior")) +
  geom_density() +
  theme_light() +
  labs(x = "Relative susceptibility of vaccinated and isolated adults", y = "Density") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey"), labels = c("1", "2", "3", "prior"))
ggsave("figures/rSAVI_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSAI), ], aes(x = rSAI, color = chain)) +
  geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSAI, na.rm = T), ],
            aes(x = x, y = d, color = "prior")) +
  geom_density() +
  theme_light() +
  labs(x = "Relative susceptibility of isolated and not vaccinated adults", y = "Density") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey"), labels = c("1", "2", "3", "prior"))
ggsave("figures/rSAI_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSC), ], aes(x = rSC, color = chain)) +
  geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSC, na.rm = T), ],
            aes(x = x, y = d, color = "prior")) +
  geom_density() +
  theme_light() +
  labs(x = "Relative susceptibility of not isolated children", y = "Density") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey"), labels = c("1", "2", "3", "prior"))
ggsave("figures/rSC_raw.pdf", width = 5, height = 4)

ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSCI), ], aes(x = rSCI, color = chain)) +
  geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSCI, na.rm = T), ],
            aes(x = x, y = d, color = "prior")) +
  geom_density() +
  theme_light() + 
  labs(x = "Relative susceptibility of isolated children", y = "Density") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey"), labels = c("1", "2", "3", "prior"))
ggsave("figures/rSCI_raw.pdf", width = 5, height = 4)


############################################
## Compare susceptibility of vaccinated + 
## isolated contacts to vaccinated not isolated 
## contacts
############################################
mcmc_posterior %>%
  filter(chain == "1") %>%
  summarise(
    `Isolated and vaccinated adults/teenagers are less susceptible than vaccinated not isolated adults/teenagers` = sum(rSAVI < rSAV)/n()
            ) %>%
  pivot_longer(everything(), names_to = "hypothesis", values_to = "Probability") %>%
  mutate(
    BayesFactor = Probability/(1-Probability),
    Evidence = log10(Probability/(1-Probability))
    ) %>%
write.table(., "tables/isolation_in_vaccinated_contacts.csv", sep=",", row.names = F)

############################################
## Parameter estimates
############################################
# Relative susceptibility 
baseline = data.frame(
  rS = "Adult",
  mean = 1,
  median = 1, 
  q025 = 1, 
  q975 = 1
)

G1 = mcmc_posterior %>%
  filter(chain == "1") %>%
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
  filter(chain == "1") %>%
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
ggsave("figures/parameter_estimates_baseline.pdf", g, width = 7.5, height = 4)


# Tables 
tab_sus = mcmc_posterior %>%
  filter(chain == "1") %>%
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
  filter(chain == "1") %>%
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
## Person-to-person transmission risk
############################################
# Absolute person-to-person risk of infection
# according to household size and 
pInf_input = mcmc_posterior %>%
  filter(chain == "1") %>%
  select(beta, delta, rSAI, rSAVI, rSAV, rSC, rSCI, rInfVac)

# Compute person-to-person transmission probability for
# vaccinated and unvaccinated cases
p_to_p_transmission <- function(df, case, contact) {
  
  df$B <- with(df, beta)
  
  if (contact != "rSA") df$B <- with(df, B*eval(parse(text=contact)))
  if (grepl("Vaccinated adult", case)) df$B <- with(df, B*rInfVac)
  
  df
}

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

# Save estimates in table
pInf %>%
  mutate(
    median = round(median * 100, 1), 
    q025 = round(q025 *100, 1), 
    q975 = round(q975*100, 1), 
    contact = factor(contact, 
                     levels = c("rSA", "rSC", "rSCI", "rSAI", "rSAV", "rSAVI"), 
                     labels = c("Adult contact", "Child contact", "Isolated child contact", "Isolated adult contact", "Vaccinated adult contact", 
                                "Isolated + vaccinated\nadult contact"))) %>%
  mutate(value = paste0(median, " (", q025, "-", q975, ")")) %>%
  pivot_wider(-c(median, q025, q975), values_from = value, names_from = case) %>%
  write.table(., "tables/transmission_risk_baseline.csv", sep=",", row.names = F)

# Figure
pInf %>%
  mutate(contact = factor(contact, 
                          levels = c("rSA", "rSC", "rSCI", "rSAI", "rSAV", "rSAVI"), 
                          labels = c("Adult contact", "Child contact", "Isolated child contact", "Isolated adult contact", "Vaccinated adult contact", 
                                     "Isolated + vaccinated\nadult contact"))) %>%
  ggplot(., aes(x = contact, y = median, ymin = q025, ymax = q975, col = case)) +
  geom_pointrange(position = position_dodge(width = 0.4), size = 0.8, fatten = 2.5) +
  theme_light() +
  theme(strip.text.x = element_text(size = 12), 
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#5B0E2D", "#FFA781")) +
  labs(y = "Probability of person-to-person transmission", x = "", col = "", title = "")
ggsave("figures/transmission_risk_baseline.pdf", width = 6.5, height = 4.5)


############################################
## Model adequation
############################################
# Household id of households with 1 index case
hhids_1_index = sapply(unique(bdd$hhid), function(x) sum(bdd$index[bdd$hhid %in% x]))
hhids_1_index = unique(bdd$hhid)[hhids_1_index == 1]

obsData = bdd %>% 
  filter(index == 0, hhid %in% hhids_1_index) %>% 
  group_by(hhsize) %>% 
  summarise(SAR = sum(infectionStatus > 0) / n()) %>%
  filter(hhsize < 12) %>%
  data.frame()

# Expected SAR
sar_exp = read.csv("results/expectedSAR_full_database_2doses.csv") %>%
  filter(hhsize < 12) %>%
  group_by(sim, hhsize) %>%
  summarise(sar = nInfContacts / nContacts) %>%
  group_by(hhsize) %>%
  summarise(
    median = median(sar), 
    q025 = quantile(sar, 0.025),
    q975 = quantile(sar, 0.975)
  )

ggplot(sar_exp[sar_exp$hhsize <= 5, ]) +
  geom_pointrange(aes(x = as.factor(hhsize), 
                      y = median, 
                      ymin = q025, 
                      ymax = q975),
                  col = "#358597") +
  geom_point(data = obsData[obsData$hhsize <= 5, ], 
             aes(x = as.factor(hhsize), y = SAR), 
             col = "black",
             shape = 15
  ) +
  theme_light() +
  ylim(c(0,1)) +
  labs(x = "Number of household members", y = "Secondary attack rate", col = "", shape = "")
ggsave(paste0("figures/model_adequacy.pdf"), height = 3.5, width = 4)

