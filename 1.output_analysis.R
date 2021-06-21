###################################################
##                  ANALYSIS
##            BASIC STATISTICS
##
## author: Maylis Layan
## creation date: 2021/04/14
## last update: 
###################################################

rm(list = ls())
# setwd("/mnt/gaia/MMMICovid/maylayan/Israel/")
setwd("V:/maylayan/Israel/")

library(tidyverse)
library(gt)
library(gridExtra)
# source("R/stat_functions.R")

####################################################
## Load and prepare data
####################################################
# Load
bdd = read.csv("Data/2021_05_14_cleaned_data.csv", 
               stringsAsFactors = F, header = T, encoding = "UTF-8")
dim(bdd)

# Households
hh = unique(bdd$hhid)

############################################
## Load outputs
############################################
burnIn = 1000

mcmc_posterior = data.frame()
mcmc_LL = data.frame()
mcmc_accept = data.frame()

parameterNames = c("alpha", "beta", "delta", "rSAI", "rSAV", "rSAVI", "rSC", "rSCI", "rInfVac")

for (vacDef in c("1dose", "2doses")) {
  for (f in list.files(paste0("Results/", vacDef, "/strict"), pattern = "^mcmc_.*_10_1.0_1.0_0.6_[0-9].txt$", full.names = T)) {
    
    # Get run info
    database = regmatches(f, gregexpr("(?<=mcmc_)[a-z]+_[a-z]+(?=_)", f, perl = TRUE))[[1]]
    model = gsub("Results.*/mcmc_[a-z]+_[a-z]+_|_[0-9]{1}.*", "", f)
    asympPeriod = gsub("Results.*/mcmc_[a-zA-Z_]+_[0-9]+_|_[0-9]{1}\\.[0-9]{1}_.*.txt", "", f)
    frequencyDependent = gsub("Results.*/mcmc_[a-zA-Z_]+_|_[0-9]{2}.*.txt", "", f)
    sdLnorm = as.numeric(gsub("Results.*/mcmc_[a-zA-Z_]+_[0-9]{1}_[0-9]{2}_|_[0-9]\\.[0-9]_[0-9]\\.[0-9]_[0-9].txt", "", f))
    chain = gsub("Results.*/.*_|.txt", "", f)
    
    # Load output data
    out_temp = read.table(f, sep = " ", header = T)
    
    # LogLik
    LL_temp = subset(out_temp, subset = iteration >= burnIn, select = c(iteration, logLik)) 
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
    
    posterior_temp = out_temp[out_temp$iteration >= burnIn, c("iteration", params[params != "data"]) ]
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

# ############################################
# ## Acceptance rates
# ############################################
# mcmc_accept = mcmc_accept %>%
#   mutate(hhSize = factor(hhSize, levels = c("0", "1"), labels = c("density dependent", "frequency dependent")),
#          parameter = factor(parameter, c(parameterNames, "data")))
# 
# ggplot(mcmc_accept[!is.na(mcmc_accept$rate) & mcmc_accept$vacDef == "1dose", ],
#        aes(x = parameter, y = rate, col = chain)) +
#   facet_grid(database ~ model) +
#   geom_jitter() +
#   ylim(c(0,1)) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
#   labs(x = "", y = "Acceptance rate")
# ggsave("Figures/1dose/test3/acceptance_rate_known_outcome.pdf", height = 5, width = 8)
# 
# ggplot(mcmc_accept[!is.na(mcmc_accept$rate) & mcmc_accept$vacDef == "2doses", ],
#        aes(x = parameter, y = rate, col = chain)) +
#   facet_grid(database ~ model) +
#   geom_jitter() +
#   ylim(c(0,1)) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
#   labs(x = "", y = "Acceptance rate")
# ggsave("Figures/2doses/test3/acceptance_rate_known_outcome.pdf", height = 5, width = 8)
# 
# mcmc_accept %>%
#   filter(!is.na(rate)) %>%
#   group_by(sdLnorm, model, database, parameter) %>%
#   summarise(mean = mean(rate)) %>%
#   data.frame()
#    
# ############################################
# ## Mixing
# ############################################
# params = params[1:9]
# 
# for (v in c("2doses", "1dose")){
#   for (d in c("full_database", "known_outcome")) {
#     pdf(paste0("Figures/", v, "/test3/mixing_", d, "_baseline.pdf"), width = 18, height = 8)
#     par(mfrow = c(3, length(params) + 1), mar = c(4,2,3,2))
#     for (r in 1:3) {
#       plot(
#         mcmc_LL$iteration[mcmc_LL$chain == r & mcmc_LL$vacDef == v & mcmc_LL$database == d], 
#         mcmc_LL$logLik[mcmc_LL$chain == r & mcmc_LL$vacDef == v & mcmc_LL$database == d],
#         type = "l", 
#         xlab = "iteration", 
#         ylab = "", 
#         main = "LogLik"
#       )
#       
#       for (p in params) {
#         plot(
#           mcmc_posterior$iteration[mcmc_posterior$chain == r & mcmc_posterior$vacDef == v & mcmc_posterior$database == d], 
#           mcmc_posterior[[p]][mcmc_posterior$chain == r & mcmc_posterior$vacDef == v & mcmc_posterior$database == d],
#           type = "l", 
#           xlab = "iteration", 
#           ylab = "", 
#           main = p
#         )
#       }
#     }
#     dev.off()
#   }
# }
# 
# 
# ############################################
# ## Verify concordance between chains
# ############################################
# logNormalPrior = data.frame(d = dlnorm(seq(0, 10, 0.01), 1), x=seq(0, 10, 0.01))
# 
# for (v in c("1dose", "2doses")) {
#   ggplot(mcmc_posterior, aes(x = alpha, color = chain)) +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v,"/test3/alpha_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior, aes(x = beta, color = chain)) +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v,"/test3/beta_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$delta), ], aes(x = delta, color = chain)) +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/delta_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$rInfVac), ], aes(x = rInfVac, color = chain)) +
#     geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rInfVac, na.rm = T), ],
#               aes(x = x, y = d), color = "grey") +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/rInfVac_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSAV), ], aes(x = rSAV, color = chain)) +
#     geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSAV, na.rm = T), ],
#               aes(x = x, y = d), color = "grey") +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/rSAV_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSAVI), ], aes(x = rSAVI, color = chain)) +
#     geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSAVI, na.rm = T), ],
#               aes(x = x, y = d), color = "grey") +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/rSAVI_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSAI), ], aes(x = rSAI, color = chain)) +
#     geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSAI, na.rm = T), ],
#               aes(x = x, y = d), color = "grey") +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/rSAI_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSC), ], aes(x = rSC, color = chain)) +
#     geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSC, na.rm = T), ],
#               aes(x = x, y = d), color = "grey") +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/rSC_raw.pdf"), width = 8, height = 5)
# 
#   ggplot(mcmc_posterior[!is.na(mcmc_posterior$rSCI), ], aes(x = rSCI, color = chain)) +
#     geom_line(data = logNormalPrior[logNormalPrior$x <= max(mcmc_posterior$rSCI, na.rm = T), ],
#               aes(x = x, y = d), color = "grey") +
#     geom_density() +
#     facet_grid(database ~ model) +
#     theme_light()
#   ggsave(paste0("Figures/", v, "/test3/rSCI_raw.pdf"), width = 8, height = 5)
# }

############################################
## Compare categories - full database - 2doses
############################################
mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "full_database", vacDef == "2doses") %>%
  summarise(child = sum(rSCI < rSC)/n(), 
            IVltV = sum(rSAVI < rSAV)/n(),
            IltV = sum(rSAI < rSAV)/n(),
            IVltI = sum(rSAVI < rSAI)/n(),
            ) %>%
  pivot_longer(everything(), names_to = "hypothesis", values_to = "p") %>%
  mutate(
    bf = p/(1-p),
    evidence = log10(p/(1-p))
    )


mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "known_outcome", vacDef == "2doses") %>%
  summarise(child = sum(rSCI < rSC)/n(), 
            IVltV = sum(rSAVI < rSAV)/n(),
            IltV = sum(rSAI < rSAV)/n(),
            IVltI = sum(rSAVI < rSAI)/n(),
  ) %>%
  pivot_longer(everything(), names_to = "hypothesis", values_to = "p") %>%
  mutate(
    bf = p/(1-p),
    evidence = log10(p/(1-p))
  )


############################################
## Fig 2 - 2 doses - Full database
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
  geom_pointrange(position = position_dodge(0.6), col= "#358597") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 10)) +
  # scale_y_log10(breaks = seq(0.1, 1.1, 0.2), minor_breaks = seq(0.1,1.1, 0.1), limits= c(0.06,1.1)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1.1)) +
  labs(y = "", x = "", title = "Relative infectivity")

g = grid.arrange(G1,G2, widths = c(2, 1))
ggsave("Paper/Fig2/fig2.pdf", g, width = 7.5, height = 4)


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
                     c("Adult", "Isolated adult", "Vaccinated adult", "Isolated + vaccinated adult", "Child", "Isolated child"))
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
  write.table(., "Paper/Supplementary_Tables/parameter_estimates_baseline.csv", sep=",", col.names = F, row.names = F)

############################################
## Fig 3 - 2 doses - full database
############################################
# Absolute person-to-person risk of infection
# according to household size and 
pInf_input = mcmc_posterior %>%
  filter(chain == "1", model == "rS_rInfVac", database == "full_database", vacDef == "2doses") %>%
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
  write.table(., "Paper/Supplementary_Tables/transmission_risk.csv", sep=",", row.names = F)


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
ggsave("Paper/Fig3/fig3_p_to_p.pdf", width = 6.5, height = 4.5)
