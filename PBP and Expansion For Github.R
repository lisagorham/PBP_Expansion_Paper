#This is the code for my PBP outcomes/expansion paper (October 2026)

#basic R setup ------------------------------------------------------------------------
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(writexl)
library(ggseg)
library(ggsegGordon)
library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(cmdstanr)
library(marginaleffects)
library(forcats)
library(emmeans)
library(modelsummary)
library(mgcv)
library(lme4)
library(lmerTest)
library(purrr)
library(ggsegGlasser)
library(RColorBrewer)
library(patchwork)
library(lme4)

setwd("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/R files")

#Load in the data --------------------------------------------------------------------------------

#grabbing outcome data for the entire broader cohort 
parent_PBP <- read.csv("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/finaldata/PBP_score_parent.csv")
self_PBP <- read.csv("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/finaldata/Child_report_PBP_scores.csv")
teacher_PBP <- read.csv("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/finaldata/Teacher_PBP_scores.csv")

#deriving the autism, anxiety, and adhd domains
parent_PBP <- parent_PBP %>%
  mutate(
    PBP_Parent = PBP_score_p,
    scared_zscore = as.numeric(scale(scared_totalanxiety_raw)),
    cbcl_zscore = as.numeric(scale(cbcl618_intprob_tscore)),
    hyperactivity_zscore = as.numeric(scale(con3_par_hyp_tscore)),
    inattention_zscore = as.numeric(scale(con3_par_inatt_tscore)),
    srs_zscore = as.numeric(scale(srs_par_total_tscore)),
    anxiety_parent_v1 = (scared_zscore + cbcl_zscore)/2,
    anxiety_parent_v2 = na_if(rowMeans(pick(scared_zscore, cbcl_zscore), na.rm = TRUE), NaN),
    adhd_parent = (hyperactivity_zscore + inattention_zscore)/2,
    autism_parent = srs_zscore,
    
  )



parent_selected <- parent_PBP %>%
  select(c(
    recruitmentid,
    PBP_Parent,
    anxiety_parent_v1,
    anxiety_parent_v2,
    adhd_parent,
    autism_parent
  ))

self_PBP <- self_PBP %>%
  mutate(
    nepsysocialz = as.numeric(scale(Znepsy_total_tscore_revs)),
    nepsytomz = as.numeric(scale(Znepsy_tomtotal_raw_revs)),
    anxiety_child = as.numeric(scale(Zbpmy_int_tscore)),
    adhd_child = as.numeric(scale(Zbpmy_att_tscore)),
    autism_child_v1 = (nepsysocialz + nepsytomz)/2,
    autism_child_v2 = na_if(rowMeans(pick(nepsysocialz, nepsytomz), na.rm = TRUE), NaN)
  )

self_selected <- self_PBP %>%
  select(c(
    recruitmentid,
    PBP_child,
    anxiety_child,
    adhd_child,
    autism_child_v1,
    autism_child_v2
  ))

teacher_PBP <- teacher_PBP %>%
  mutate(
    hyperactivity_zscoreteach = as.numeric(scale(con3_tea_hyp_tscore)),
    inattention_zscoreteach = as.numeric(scale(con3_tea_inatt_tscore)),
    anxiety_teacher = as.numeric(scale(trf618_intprob_tscore)),
    autism_teacher = as.numeric(scale(srs_tea_total_tscore)),
    adhd_teacher_v1 = (hyperactivity_zscoreteach + inattention_zscoreteach)/2,
    adhd_teacher_v2 = na_if(rowMeans(pick(hyperactivity_zscoreteach, inattention_zscoreteach), na.rm = TRUE), NaN)
  )

teacher_selected <- teacher_PBP %>%
  select(c(
    recruitmentid,
    PBP_teacher,
    anxiety_teacher,
    autism_teacher,
    adhd_teacher_v1,
    adhd_teacher_v2
  ))

outcomes <- parent_selected %>% 
  full_join(self_selected) %>%
  full_join(teacher_selected)

#grabbing demographic information for the 121 subjects included in analyses
personleveldata <- read_xlsx("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/personleveldata_10.20.26_revised.xlsx")
personleveldata <- personleveldata %>% left_join(outcomes)
wholebrainexp <- read_xlsx("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/wholebrainexpnumbers.xlsx")
personleveldata <- personleveldata %>% left_join(wholebrainexp)
personleveldata$wholebrainexp_zscore <- as.numeric(scale(personleveldata$whole_brain_expansion_avg))
personleveldata$PBP_Parent_zscore <- as.numeric(scale(personleveldata$PBP_Parent))
personleveldata$child_PBP_zscore <- as.numeric(scale(personleveldata$PBP_child))
personleveldata$teacher_PBP_zscore <- as.numeric(scale(personleveldata$PBP_teacher))
personleveldata$anxiety_parent_zscore <- as.numeric(scale(personleveldata$anxiety_parent_v2))
personleveldata$adhd_parent_zscore <- as.numeric(scale(personleveldata$adhd_parent))
personleveldata$autism_parent_zscore <- as.numeric(scale(personleveldata$autism_parent))
personleveldata$anxiety_child_zscore <- as.numeric(scale(personleveldata$anxiety_child))
personleveldata$adhd_child_zscore <- as.numeric(scale(personleveldata$adhd_child))
personleveldata$autism_child_zscore <- as.numeric(scale(personleveldata$autism_child_v2))

#load in the glasser expansiondata
glasser <- read.csv("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/updated data from aidan/glassernetworks.csv")
colnames(glasser)[1] = "ID_910_USE"
long_glasser <- glasser %>%
  pivot_longer(
    cols = -ID_910_USE,
    names_to = "Network",
    values_to = "Network_Expansion"
  )

completeglasser <- personleveldata %>% left_join(long_glasser)
completeglasser$Network_Expansion_Zscore <- as.numeric(scale(completeglasser$Network_Expansion))

#Figure 1 code (group averages and plotting networks on brain) ---------------------------------------------------------------
#first, create the lookup table
glasser_atlas <- ggsegGlasser::glasser
glasser_atlas

lookuptable <- read_xlsx("C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/glasser_network_lookup.xlsx")
lookuptable$NETWORK <- gsub("-", ".", lookuptable$NETWORK)
lookuptable$NETWORK <- gsub(
  "Dorsal.Attention",
  "Dorsal.attention",
  lookuptable$NETWORK
)
lookuptable <- lookuptable %>%
  separate(
    GLASSERLABELNAME,
    into = c("hemi", "region", "remove"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    hemi = ifelse(hemi == "L", "left", "right")
  ) %>%
  select(-remove)
lookuptable <- lookuptable %>% select(c(NETWORK, hemi, region))
colnames(lookuptable)[1] = "Network"

#Make the 3 datasets for plotting
FT <- completeglasser %>% filter(group == "FT")
VPT_injury <- completeglasser %>% filter(group == "VPT_injury")
VPT_noinjury <- completeglasser %>% filter(group == "VPT_noinjury")

FT_plotting <- lookuptable %>% left_join(FT)
VPT_injury_plotting <- lookuptable %>% left_join(VPT_injury)
VPT_noinjury_plotting <- lookuptable %>% left_join(VPT_noinjury)


#use the completeglasser$Network_Expansion variable



#plot the FT kids average expansion
ggseg(FT_plotting,
      colour = "black", #changes the colors of the outlines of the ROIs
      atlas = glasser_atlas, 
      position = "dispersed",
      mapping = aes(fill = Network_Expansion)) +
  theme_classic() +
  scale_fill_gradientn(
    colours = c(
      "#3B4CC0",  # blue
      "#688BEA",  # periwinkle
      "#2FB7A5",  # teal
      "#7AD151",  # green
      "#FDE725",  # yellow
      "#F89C20",  # orange
      "#E64B35",  # red
      "#B12A90"   # magenta
    ),
    limits = c(0.94, 2.05)
  ) +
  labs(title = "Mean Expansion per Network for FT Kids", fill = "Expansion")

ggseg(VPT_noinjury_plotting,
      colour = "black", #changes the colors of the outlines of the ROIs
      atlas = glasser_atlas, 
      position = "dispersed",
      mapping = aes(fill = Network_Expansion)) +
  theme_classic() +
  scale_fill_gradientn(
    colours = c(
      "#3B4CC0",  # blue
      "#688BEA",  # periwinkle
      "#2FB7A5",  # teal
      "#7AD151",  # green
      "#FDE725",  # yellow
      "#F89C20",  # orange
      "#E64B35",  # red
      "#B12A90"   # magenta
    ),
    limits = c(0.94, 2.05)
  ) +
  labs(title = "Mean Expansion per Network for Non-Injured VPT Kids", fill = "Expansion")

ggseg(VPT_injury_plotting,
      colour = "black", #changes the colors of the outlines of the ROIs
      atlas = glasser_atlas, 
      position = "dispersed",
      mapping = aes(fill = Network_Expansion)) +
  theme_classic() +
  scale_fill_gradientn(
    colours = c(
      "#3B4CC0",  # blue
      "#688BEA",  # periwinkle
      "#2FB7A5",  # teal
      "#7AD151",  # green
      "#FDE725",  # yellow
      "#F89C20",  # orange
      "#E64B35",  # red
      "#B12A90"   # magenta
    ),
    limits = c(0.94, 2.05)
  ) +
  labs(title = "Mean Expansion per Network for Injured VPT Kids", fill = "Expansion")+
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))

#showing the networks on the brain
ggseg(FT_plotting,
      colour = "black",
      atlas = glasser_atlas, 
      position = "dispersed",
      mapping = aes(fill = Network)) +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Visual1" = "#E41A1C",
      "Visual2" = "#FF7F00",
      "Somatomotor" = "#FFD92F",
      "Cingulo.Opercular" = "#4DAF4A",
      "Dorsal.attention" = "#00A6A6",
      "Language" = "#377EB8",
      "Frontoparietal" = "#6A3D9A",
      "Auditory" = "#984EA3",
      "Default" = "#F781BF",
      "Posterior.Multimodal" = "#A65628",
      "Ventral.Multimodal" = "#66C2A5",
      "Orbito.Affective" = "#1F78B4"
    )
  ) +
  labs(
    title = "Glasser Networks",
    fill = "Network"
  )

#Demographic analyses ---------------------------------------------------------------------------------------------

personleveldata %>%
  group_by(group) %>%
  summarise(
    regvis_scanpma = paste0(round(mean(regvis_scanpma, na.rm = TRUE), 1),
                            " (", round(sd(regvis_scanpma, na.rm = TRUE), 1), ")"),
    poverty = paste0(round(mean(neighborhood_natlcentile_revised, na.rm = TRUE), 1),
                     " (", round(sd(neighborhood_natlcentile_revised, na.rm = TRUE), 1), ")"),
    age910 = paste0(round(mean(mrinotes_age_yr, na.rm = TRUE), 1),
                    " (", round(sd(mrinotes_age_yr, na.rm = TRUE), 1), ")"),
    Parent_PBP = paste0(round(mean(PBP_Parent, na.rm = TRUE), 1),
                 " (", round(sd(PBP_Parent, na.rm = TRUE), 1), ")"),
    Child_PBP = paste0(round(mean(PBP_child, na.rm = TRUE), 1),
                        " (", round(sd(PBP_child, na.rm = TRUE), 1), ")"),
    Teacher_PBP = paste0(round(mean(PBP_teacher, na.rm = TRUE), 1),
                        " (", round(sd(PBP_teacher, na.rm = TRUE), 1), ")"),
    BirthSA = paste0(round(mean(birth_totalSA, na.rm = TRUE), 1),
                     " (", round(sd(birth_totalSA, na.rm = TRUE), 1), ")"),
    age910SA = paste0(round(mean(age910_totalSA, na.rm = TRUE), 1),
                      " (", round(sd(age910_totalSA, na.rm = TRUE), 1), ")")
  )

personleveldata %>%
  group_by(group) %>%
  summarise(
    gestationalage = paste0(round(mean(birth_ga, na.rm = TRUE), 1), ",",
                            min(birth_ga, na.rm = TRUE), ",",
                            max(birth_ga, na.rm = TRUE),
                            " (", round(sd(birth_ga, na.rm = TRUE), 1), ")") )




#for continuous variables
summary(aov(age910_totalSA ~ group, data = personleveldata))

#for categorical variables
chisq.test(table(personleveldata$sex, personleveldata$group))

table(personleveldata$sex, personleveldata$group)

#checking to see if sex is related to expansion or PBP
t.test(whole_brain_expansion_avg ~ sex, data = personleveldata)

#Surface area at endpoints models --------------------------------------------------------

birthSAnoint <- brm(wholebrainbirthSA_zscore ~ PBP_Parent_zscore + group + poverty_zscore + pma_zscore + sex,
                    data = personleveldata,
                    family = gaussian(),
                    chains = 4, iter = 4000, warmup = 2000, cores = 4,
                    backend = "cmdstanr",
                    control = list(
                      adapt_delta = 0.99,
                      max_treedepth = 15
                    ), file = "birthSA_noint_10.21.26"
)

summary(birthSAnoint)

personleveldata$age910_zscore <- as.numeric(scale(personleveldata$mrinotes_age_yr))

age910SAnoint <- brm(wholebrainage910SA_zscore ~ PBP_Parent_zscore + group + poverty_zscore + age910_zscore + sex,
                     data = personleveldata,
                     family = gaussian(),
                     chains = 4, iter = 4000, warmup = 2000, cores = 4,
                     backend = "cmdstanr",
                     control = list(
                       adapt_delta = 0.99,
                       max_treedepth = 15
                     ), file = "age910SA_noint_10.21.26"
)

summary(age910SAnoint)

age910SAint <- brm(wholebrainage910SA_zscore ~ PBP_Parent_zscore*group + poverty_zscore + age910_zscore + sex,
                   data = personleveldata,
                   family = gaussian(),
                   chains = 4, iter = 4000, warmup = 2000, cores = 4,
                   backend = "cmdstanr",
                   control = list(
                     adapt_delta = 0.99,
                     max_treedepth = 15
                   ), file = "age910SA_int_10.21.26"
)

summary(age910SAint)





#question 1- is there a whole brain relationship with symptoms? --------------------------------------------------------------

wholebrainnoint <- lm(wholebrainexp_zscore ~ PBP_Parent_zscore + group + poverty_zscore + pma_zscore, data = personleveldata)
summary(wholebrainnoint)

#bayesian version
wholebrainnoint_bayes <- brm(wholebrainexp_zscore ~ PBP_Parent_zscore + group + poverty_zscore + pma_zscore,
                             data = personleveldata,
                             family = gaussian(),
                             chains = 4, iter = 4000, warmup = 2000, cores = 4,
                             backend = "cmdstanr",
                             control = list(
                               adapt_delta = 0.99,
                               max_treedepth = 15
                             ), file = "wholebrainmodelnoint_bayes_10.20.26"
)

summary(wholebrainnoint_bayes)

#trying with an interaction to see if it compares between groups
wholebrainint <- lm(wholebrainexp_zscore ~ PBP_Parent_zscore*group + poverty_zscore + pma_zscore, data = personleveldata)
summary(wholebrainint)

#bayesian version
wholebrainint_bayes <- brm(wholebrainexp_zscore ~ PBP_Parent_zscore*group + poverty_zscore + pma_zscore,
                             data = personleveldata,
                             family = gaussian(),
                             chains = 4, iter = 4000, warmup = 2000, cores = 4,
                             backend = "cmdstanr",
                             control = list(
                               adapt_delta = 0.99,
                               max_treedepth = 15
                             ), file = "wholebrainmodelint_bayes_10.20.26"
)

summary(wholebrainint_bayes)

#plotting raw data
p_raw <- ggplot(
  personleveldata,
  aes(
    x = PBP_Parent_zscore,
    y = 2^whole_brain_expansion_avg,
    fill = group
  )
) +
  geom_point(
    shape = 21,
    alpha = 0.7,
    size = 2.5
  ) +
  labs(
    x = "PBP (z-score)",
    y = "Whole brain expansion (A2/A1)",
    fill = "Group"
  ) +
  theme_classic()+
  scale_x_continuous(
    limits = c(-2, 4),
    breaks = seq(-2, 4, by = 1)
  )+
  scale_y_continuous(
    limits = c(2.1, 3.8),
    breaks = seq(2.1, 3.8, by = 0.2)
  )
p_raw

#plotting predicted effects for the three groups
mean_wholebrain <- mean(personleveldata$whole_brain_expansion_avg, na.rm = TRUE)
sd_wholebrain <- sd(personleveldata$whole_brain_expansion_avg, na.rm = TRUE)

ce <- conditional_effects(
  wholebrainnoint_bayes,
  effects = "PBP_Parent_zscore:group"
)

ce_data <- ce$`PBP_Parent_zscore:group`
ce_data <- ce_data %>%
  mutate(
    estimate_raw = (estimate__ * sd_wholebrain) + mean_wholebrain,
    lower_raw = (lower__ * sd_wholebrain) + mean_wholebrain,
    upper_raw = (upper__ * sd_wholebrain) + mean_wholebrain,
    
    estimate_ratio = 2^estimate_raw,
    lower_ratio = 2^lower_raw,
    upper_ratio = 2^upper_raw
  )
p_model <- ggplot(
  ce_data,
  aes(
    x = PBP_Parent_zscore,
    y = estimate_ratio,
    color = group,
    fill = group
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_ratio,
      ymax = upper_ratio
    ),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  labs(
    x = "PBP (z-score)",
    y = "Predicted whole brain expansion (A2/A1)",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic()+
  scale_x_continuous(
    limits = c(-2, 4),
    breaks = seq(-2, 4, by = 1)
  )+
  scale_y_continuous(
    limits = c(2.1, 3.8),
    breaks = seq(2.1, 3.8, by = 0.2)
  )

p_model

#showing both (old version)
p_raw + p_model

#generating predictions
pred <- as.data.frame(fitted(
  wholebrainnoint_bayes,
  newdata = data.frame(
    PBP_Parent_zscore = c(-1.5, 3),
    group = "VPT_noinjury",
    poverty_zscore = 0,
    pma_zscore = 0
  ),
  summary = TRUE
))

pred

pred$estimate_raw <- pred$Estimate * sd_wholebrain + mean_wholebrain
pred$lower_raw <- pred$Q2.5 * sd_wholebrain + mean_wholebrain
pred$upper_raw <- pred$Q97.5 * sd_wholebrain + mean_wholebrain

pred$estimate_ratio <- 2^pred$estimate_raw
pred$lower_ratio <- 2^pred$lower_raw
pred$upper_ratio <- 2^pred$upper_raw

pred



#determining what a z score of 0 means for poverty and PMA
mean(personleveldata$neighborhood_natlcentile_revised)
mean(personleveldata$regvis_scanpma)


#question 2- does the effect size differ by network? --------------------------------------------

question2noint <- brm(Network_Expansion_Zscore ~ PBP_Parent_zscore + group + poverty_zscore + pma_zscore +
                            (1 | recruitmentid) + (PBP_Parent_zscore | Network),
                          data = completeglasser,
                          family = gaussian(),
                          chains = 4, iter = 4000, warmup = 2000, cores = 4,
                          backend = "cmdstanr",
                          control = list(
                            adapt_delta = 0.99,
                            max_treedepth = 15
                          ), file = "question2noint_10.20.26"
)

summary(question2noint)

#pulling out effects per network
question2noint_posterior <- question2noint %>%
  spread_draws(
    b_PBP_Parent_zscore,
    r_Network[Network, term]
  ) %>%
  filter(term == "PBP_Parent_zscore") %>%
  mutate(
    full_slope = b_PBP_Parent_zscore + r_Network
  )

question2noint_draws <- question2noint_posterior %>%
  group_by(Network) %>%
  mutate(median_slope = median(full_slope)) %>%
  ungroup() %>%
  mutate(Network = forcats::fct_reorder(Network, median_slope))

question2noint_summary <- question2noint_posterior %>%
  group_by(Network) %>%
  median_qi(full_slope, .width = 0.95) %>%
  ungroup() %>% 
  mutate(
    crosses_zero = .lower <= 0 & .upper >= 0
  )

question2noint_draws <- question2noint_draws %>%
  left_join(question2noint_summary %>% select(Network, crosses_zero),
            by = "Network")

#plotting the effect size error bars
ggplot(question2noint_summary,
       aes(x = full_slope,
           y = reorder(Network, full_slope))) +
  geom_errorbar(
    aes(xmin = .lower,
        xmax = .upper),
    orientation = "y",
    height = 0.2
  ) +
  geom_point(
    aes(fill = crosses_zero),
    shape = 21,
    color = "black",
    size = 6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "red",
    linewidth = 1.2
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "white"
    )
  ) +
  labs(
    x = "Effect Size for Relationship Between Expansion and PBP Symptoms (β)\n(A strongly negative effect means that greater expansion is associated with less impairing symptoms)",
    y = "Network"
  ) +
  theme_classic()+
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = seq(-0.3, 0.1, by = 0.02)
  )

#Doing pairwise comparisons

#first make the posterior wide
wideposterior <- question2noint_draws %>%
  select(.draw, Network, full_slope) %>%
  pivot_wider(names_from = Network,
              values_from = full_slope)

#make a list of all of the pairs of networks to compare 
networks <- setdiff(names(wideposterior), ".draw")
pairs <- combn(networks, 2, simplify = FALSE)

#calculate posterior difference for each pair of networks 
pairwise_results <- map_dfr(pairs, function(x){
  A <- x[1]
  B <- x[2]
  diff <- wideposterior[[A]] - wideposterior[[B]]
  tibble(
    network1 = A,
    network2 = B,
    median_difference = median(diff),
    lower95 = quantile(diff, 0.025),
    upper95 = quantile(diff, 0.975)
  )
})

#determine which differences have 95% CIs that exclude zero
significantpairs <- pairwise_results %>%
  filter((upper95 < 0 & lower95 < 0) | (upper95 > 0 & lower95 > 0))

# Order networks by posterior median effect size
network_order <- wideposterior %>%
  select(all_of(networks)) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "network",
    values_to = "median_effect"
  ) %>%
  arrange(median_effect) %>%
  pull(network)

# create matrix of posterior median differences between networks
diff_matrix <- matrix(
  NA,
  nrow = length(networks),
  ncol = length(networks),
  dimnames = list(networks, networks)
)

for (i in networks) {
  for (j in networks) {
    diff_matrix[i, j] <- median(wideposterior[[i]] - wideposterior[[j]])
  }
}

# Convert to dataframe for ggplot
diff_df <- as.data.frame(diff_matrix) %>%
  tibble::rownames_to_column("network1") %>%
  pivot_longer(
    cols = -network1,
    names_to = "network2",
    values_to = "difference"
  ) %>%
  mutate(
    network1 = factor(network1, levels = network_order),
    network2 = factor(network2, levels = network_order)
  ) %>%
  filter( #this line keeps just one triangle rather than the whole square
    as.numeric(network1) < as.numeric(network2)
  ) %>%
  mutate(
    plot_row = factor(network2, levels = network_order),
    plot_col = factor(network1, levels = network_order)
  )

#figure out where to put * for pair differences that have a 95% CI that excludes zero
stars_df <- significantpairs %>%
  mutate(
    network1 = factor(network1, levels = network_order),
    network2 = factor(network2, levels = network_order)
  ) %>%
  mutate(
    plot_col = if_else(
      as.numeric(network1) < as.numeric(network2),
      as.character(network1),
      as.character(network2)
    ),
    plot_row = if_else(
      as.numeric(network1) < as.numeric(network2),
      as.character(network2),
      as.character(network1)
    ),
    label = "*"
  ) %>%
  mutate(
    plot_col = factor(plot_col, levels = network_order),
    plot_row = factor(plot_row, levels = network_order)
  )


ggplot(
  diff_df,
  aes(
    x = plot_col,
    y = plot_row,
    fill = difference
  )
) +
  geom_tile(
    color = "gray70",
    linewidth = 0.5
  ) +
  geom_text(
    data = stars_df,
    aes(
      x = plot_col,
      y = plot_row,
      label = label
    ),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    name = "Column - Row\n(effect size)"
  ) +
  scale_x_discrete(limits = network_order) +
  scale_y_discrete(limits = rev(network_order)) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + labs(x = "", y = "")


#making the results table for posthoc comparisons
#making table of differences for results
pairwise_table <- pairwise_results %>%
  mutate(
    cell = sprintf(
      "%.3f [%.3f, %.3f]",
      median_difference,
      lower95,
      upper95
    )
  )


results_matrix <- matrix(
  "",
  nrow = length(network_order),
  ncol = length(network_order),
  dimnames = list(network_order, network_order)
)

for (i in seq_along(network_order)) {
  for (j in seq_along(network_order)) {
    
    # Diagonal
    if (i == j) {
      results_matrix[i, j] <- "—"
      
      # LOWER triangle only
    } else if (i > j) {
      
      row_network <- network_order[i]
      col_network <- network_order[j]
      
      result <- pairwise_results %>%
        filter(
          (network1 == row_network & network2 == col_network) |
            (network1 == col_network & network2 == row_network)
        )
      
      if (nrow(result) == 1) {
        
        difference <- result$median_difference
        lower <- result$lower95
        upper <- result$upper95
        
        # Convert to ROW - COLUMN
        if (result$network1 != row_network) {
          difference <- -difference
          lower <- -result$upper95
          upper <- -result$lower95
        }
        
        results_matrix[i, j] <- sprintf(
          "%.3f [%.3f, %.3f]",
          difference,
          lower,
          upper
        )
      }
    }
  }
}

results_table <- as.data.frame(results_matrix) %>%
  tibble::rownames_to_column("Network")

write.csv(results_table, "C:/Users/lisa.gorham/Box/2023 PhD Files/Aim 3 Outcomes/pairwiseresults.csv")

#plotting median effects on brain for panel C
datasetforplotting <- lookuptable %>% left_join(question2noint_summary)

ggseg(datasetforplotting,
      colour = "black", #changes the colors of the outlines of the ROIs
      atlas = glasser_atlas, 
      position = "stacked",
      mapping = aes(fill = ifelse(.upper > 0, NA, full_slope))) +
  theme_classic() +
  scale_fill_gradient(
    low = "red",
    high = "lightgrey",
    na.value = "white"
  ) +
  labs(title = "Median Effect of Expansion in Each Glasser Network on PBP Symptoms", fill = "Effect Size")

#trying the model with an interaction to check
question2int <- brm(Network_Expansion_Zscore ~ PBP_Parent_zscore*group + poverty_zscore + pma_zscore +
                        (1 | recruitmentid) + (PBP_Parent_zscore | Network),
                      data = completeglasser,
                      family = gaussian(),
                      chains = 4, iter = 4000, warmup = 2000, cores = 4,
                      backend = "cmdstanr",
                      control = list(
                        adapt_delta = 0.99,
                        max_treedepth = 15
                      ), file = "question2int_10.20.26"
)

summary(question2int)
#neither interaction is significant










#Sensitivity analyses- self and teacher report data ---------------------------------------------

cor.test(personleveldata$PBP_child, personleveldata$PBP_teacher)


selfnoint <- brm(wholebrainexp_zscore ~ child_PBP_zscore + group + poverty_zscore + pma_zscore,
                             data = personleveldata,
                             family = gaussian(),
                             chains = 4, iter = 4000, warmup = 2000, cores = 4,
                             backend = "cmdstanr",
                             control = list(
                               adapt_delta = 0.99,
                               max_treedepth = 15
                             ), file = "selfnoint_10.21.26"
)

summary(selfnoint)


selfint <- brm(wholebrainexp_zscore ~ child_PBP_zscore*group + poverty_zscore + pma_zscore,
                 data = personleveldata,
                 family = gaussian(),
                 chains = 4, iter = 4000, warmup = 2000, cores = 4,
                 backend = "cmdstanr",
                 control = list(
                   adapt_delta = 0.99,
                   max_treedepth = 15
                 ), file = "selfint_10.21.26"
)

summary(selfint)

teachernoint <- brm(wholebrainexp_zscore ~ teacher_PBP_zscore + group + poverty_zscore + pma_zscore,
                 data = personleveldata,
                 family = gaussian(),
                 chains = 4, iter = 4000, warmup = 2000, cores = 4,
                 backend = "cmdstanr",
                 control = list(
                   adapt_delta = 0.99,
                   max_treedepth = 15
                 ), file = "teachernoint_10.21.26"
)

summary(teachernoint)

#Sensitivity analyses- anxiety adhd and autism domains ------------------------------------

anxietyparentnoint <- brm(wholebrainexp_zscore ~ anxiety_parent_zscore + group + poverty_zscore + pma_zscore,
                 data = personleveldata,
                 family = gaussian(),
                 chains = 4, iter = 4000, warmup = 2000, cores = 4,
                 backend = "cmdstanr",
                 control = list(
                   adapt_delta = 0.99,
                   max_treedepth = 15
                 ), file = "anxietyparentnoint_10.21.26"
)

summary(anxietyparentnoint)

adhdparentnoint <- brm(wholebrainexp_zscore ~ adhd_parent_zscore + group + poverty_zscore + pma_zscore,
                          data = personleveldata,
                          family = gaussian(),
                          chains = 4, iter = 4000, warmup = 2000, cores = 4,
                          backend = "cmdstanr",
                          control = list(
                            adapt_delta = 0.99,
                            max_treedepth = 15
                          ), file = "adhdparentnoint_10.21.26"
)

summary(adhdparentnoint)

autismparentnoint <- brm(wholebrainexp_zscore ~ autism_parent_zscore + group + poverty_zscore + pma_zscore,
                       data = personleveldata,
                       family = gaussian(),
                       chains = 4, iter = 4000, warmup = 2000, cores = 4,
                       backend = "cmdstanr",
                       control = list(
                         adapt_delta = 0.99,
                         max_treedepth = 15
                       ), file = "autismparentnoint_10.21.26"
)

summary(autismparentnoint)


anxietyparentint <- brm(wholebrainexp_zscore ~ anxiety_parent_zscore*group + poverty_zscore + pma_zscore,
                          data = personleveldata,
                          family = gaussian(),
                          chains = 4, iter = 4000, warmup = 2000, cores = 4,
                          backend = "cmdstanr",
                          control = list(
                            adapt_delta = 0.99,
                            max_treedepth = 15
                          ), file = "anxietyparentint_10.21.26"
)

summary(anxietyparentint)

adhdparentint <- brm(wholebrainexp_zscore ~ adhd_parent_zscore*group + poverty_zscore + pma_zscore,
                       data = personleveldata,
                       family = gaussian(),
                       chains = 4, iter = 4000, warmup = 2000, cores = 4,
                       backend = "cmdstanr",
                       control = list(
                         adapt_delta = 0.99,
                         max_treedepth = 15
                       ), file = "adhdparentint_10.21.26"
)

summary(adhdparentint)

autismparentint <- brm(wholebrainexp_zscore ~ autism_parent_zscore*group + poverty_zscore + pma_zscore,
                         data = personleveldata,
                         family = gaussian(),
                         chains = 4, iter = 4000, warmup = 2000, cores = 4,
                         backend = "cmdstanr",
                         control = list(
                           adapt_delta = 0.99,
                           max_treedepth = 15
                         ), file = "autismparentint_10.21.26"
)

summary(autismparentint)

anxietychildnoint <- brm(wholebrainexp_zscore ~ anxiety_child_zscore + group + poverty_zscore + pma_zscore,
                          data = personleveldata,
                          family = gaussian(),
                          chains = 4, iter = 4000, warmup = 2000, cores = 4,
                          backend = "cmdstanr",
                          control = list(
                            adapt_delta = 0.99,
                            max_treedepth = 15
                          ), file = "anxietychildnoint_10.21.26"
)

summary(anxietychildnoint)

adhdchildnoint <- brm(wholebrainexp_zscore ~ adhd_child_zscore + group + poverty_zscore + pma_zscore,
                         data = personleveldata,
                         family = gaussian(),
                         chains = 4, iter = 4000, warmup = 2000, cores = 4,
                         backend = "cmdstanr",
                         control = list(
                           adapt_delta = 0.99,
                           max_treedepth = 15
                         ), file = "adhdchildnoint_10.21.26"
)

summary(adhdchildnoint)

autismchildnoint <- brm(wholebrainexp_zscore ~ autism_child_zscore + group + poverty_zscore + pma_zscore,
                      data = personleveldata,
                      family = gaussian(),
                      chains = 4, iter = 4000, warmup = 2000, cores = 4,
                      backend = "cmdstanr",
                      control = list(
                        adapt_delta = 0.99,
                        max_treedepth = 15
                      ), file = "autismchildnoint_10.21.26"
)

summary(autismchildnoint)

anxietychildint <- brm(wholebrainexp_zscore ~ anxiety_child_zscore*group + poverty_zscore + pma_zscore,
                         data = personleveldata,
                         family = gaussian(),
                         chains = 4, iter = 4000, warmup = 2000, cores = 4,
                         backend = "cmdstanr",
                         control = list(
                           adapt_delta = 0.99,
                           max_treedepth = 15
                         ), file = "anxietychildint_10.21.26"
)

summary(anxietychildint)

adhdchildint <- brm(wholebrainexp_zscore ~ adhd_child_zscore*group + poverty_zscore + pma_zscore,
                      data = personleveldata,
                      family = gaussian(),
                      chains = 4, iter = 4000, warmup = 2000, cores = 4,
                      backend = "cmdstanr",
                      control = list(
                        adapt_delta = 0.99,
                        max_treedepth = 15
                      ), file = "adhdchildint_10.21.26"
)

summary(adhdchildint)

autismchildint <- brm(wholebrainexp_zscore ~ autism_child_zscore*group + poverty_zscore + pma_zscore,
                        data = personleveldata,
                        family = gaussian(),
                        chains = 4, iter = 4000, warmup = 2000, cores = 4,
                        backend = "cmdstanr",
                        control = list(
                          adapt_delta = 0.99,
                          max_treedepth = 15
                        ), file = "autismchildint_10.21.26"
)

summary(autismchildint)
