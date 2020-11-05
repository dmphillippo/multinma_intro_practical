## Example: Plaque Psoriasis

library(multinma)
library(dplyr)      # dplyr and tidyr for data manipulation
library(tidyr)
library(ggplot2)    # ggplot2 for plotting covariate distributions

options(mc.cores = parallel::detectCores())



# Load data ---------------------------------------------------------------

pso_ipd <- filter(plaque_psoriasis_ipd,
                  studyc %in% c("UNCOVER-1", "UNCOVER-2", "UNCOVER-3"))

pso_agd <- filter(plaque_psoriasis_agd,
                  studyc == "FIXTURE")

head(pso_ipd)
head(pso_agd)



# Variable transformations ------------------------------------------------

pso_ipd <- pso_ipd %>%
  mutate(# Variable transformations
         bsa = bsa / 100,
         prevsys = as.numeric(prevsys),
         psa = as.numeric(psa),
         weight = weight / 10,
         durnpso = durnpso / 10,
         # Treatment classes
         trtclass = case_when(trtn == 1 ~ "Placebo",
                              trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
                              trtn == 4 ~ "TNFa blocker"),
         # Check complete cases for covariates of interest
         complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  )

pso_agd <- pso_agd %>%
  mutate(
    # Variable transformations
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    prevsys = prevsys / 100,
    psa = psa / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                              trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
                              trtn == 4 ~ "TNFa blocker")
  )



# Complete cases ----------------------------------------------------------

sum(!pso_ipd$complete)
mean(!pso_ipd$complete)

pso_ipd <- filter(pso_ipd, complete)



# Make network ------------------------------------------------------------

## Q: Fill in the blanks to set up the network

pso_net_ipd <- set_ipd(pso_ipd,
                       study = <???>,
                       trt = <???>,
                       r = <???>,
                       trt_class = trtclass)

pso_net_agd <- set_agd_arm(pso_agd,
                           <???>
                           trt_class = trtclass)

pso_net <- combine_network(pso_net_ipd, pso_net_agd)
pso_net

# Network plot
plot(pso_net)

# We can customise the plot
plot(pso_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")



# Examine covariates for numerical integration ----------------------------

# We look at continuous covariates only

ipd_summary <- pso_ipd %>%
  group_by(studyc) %>%
  # Get mean and sd of covariates in each study
  summarise_at(vars(weight, durnpso, bsa), list(mean = mean, sd = sd, min = min, max = max)) %>%
  pivot_longer(weight_mean:bsa_max, names_sep = "_", names_to = c("covariate", ".value")) %>%
  # Assign distributions
  mutate(dist = recode(covariate,
                       bsa = "dlogitnorm",
                       durnpso = "dgamma",
                       weight = "dgamma")) %>%
  # Compute density curves
  group_by(studyc, covariate) %>%
  mutate(value = if_else(dist == "dlogitnorm",
                         list(seq(0, 1, length.out = 101)),
                         list(seq(min*0.8, max*1.2, length.out = 101)))) %>%
  unnest(cols = value) %>%
  mutate(dens = do.call(first(dist), args = list(x = value, mean = first(mean), sd = first(sd))))

# Plot histograms and assumed densities
pso_ipd %>%
  pivot_longer(c(weight, durnpso, bsa), names_to = "covariate", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = stat(density)),
                 binwidth = function(x) diff(range(x)) / nclass.Sturges(x),
                 boundary = 0,
                 fill = "grey50") +
  geom_line(aes(y = dens), data = ipd_summary,
            colour = "darkred", size = 0.5) +
  facet_wrap(~studyc + covariate, scales = "free", ncol = 3) +
  theme_multinma()



# Add numerical integration points to network -----------------------------

## Q: Fill in the blanks to add numerical integration points to the network.
##    Use the following distributions:
##      - durnpso = Gamma
##      - prevsys = Bernoulli
##      - bsa = logit-Normal
##      - weight = Gamma
##      - psa = Bernoulli

pso_net <- add_integration(pso_net,
  durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
  prevsys = distr(qbern, prob = prevsys),
  bsa = distr(qlogitnorm, <???>),
  weight = distr(<???>),
  psa = <???>,
  n_int = 1000
)

pso_net



# Fit ML-NMR model --------------------------------------------------------

## Q: Modify the code below to fit a FE ML-NMR model with a probit link, and
##    prognostic (main) effects and effect modifying interactions for all five
##    covariates: durnpso, prevsys, bsa, weight, psa

summary(normal(scale = 10))

pso_fit_FE <- nma(pso_net,
                  trt_effects = "fixed",
                  link = <???>,
                  likelihood = "bernoulli2",
                  regression = ~durnpso*.trt,
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10),
                  init_r = 0.1,
                  QR = TRUE)

print(pso_fit_FE)

# Compare prior and posterior
plot_prior_posterior(pso_fit_FE)

# Check numerical integration error
plot_integration_error(pso_fit_FE)

# Plot posterior of regression parameters
plot(pso_fit_FE,
     pars = "beta",
     stat = "halfeye",
     ref_line = 0)


# Fit RE model ------------------------------------------------------------
#
# # The following code runs  RE model, to check for residual heterogeneity
# # (possibly due to unobserved effect modifiers). We won't run this now, as
# # it takes a few minutes to run.
#
# summary(normal(scale = 10))
# summary(half_normal(scale = 2.5))
#
# pso_fit_RE <- nma(pso_net,
#                   trt_effects = "random",
#                   link = "probit",
#                   likelihood = "bernoulli2",
#                   regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
#                   class_interactions = "common",
#                   prior_intercept = normal(scale = 10),
#                   prior_trt = normal(scale = 10),
#                   prior_reg = normal(scale = 10),
#                   prior_het = half_normal(scale = 2.5),
#                   init_r = 0.1,
#                   QR = TRUE)
#
#
# print(pso_fit_RE)
#
#
# Examine divergent transitions
# pairs(pso_fit_RE, pars = c("delta[UNCOVER-2: ETN]", "d[ETN]", "tau", "lp__"))
#
# Compare prior and posterior
# plot_prior_posterior(pso_fit_RE, prior = "het")
#
# Check integration error
# plot_integration_error(pso_fit_RE)


# Compare DIC -------------------------------------------------------------

pso_dic_FE <- dic(pso_fit_FE)
pso_dic_FE
# pso_dic_RE <- dic(pso_fit_RE)
# pso_dic_RE



# Results for each study population in the network ------------------------

## Q: Produce and plot population-average relative effects and event
##    probabilities for each study in the network
##
##    (hint: use the relative_effects() and predict() functions)

# Relative effects

pso_releff_FE <- <???>
pso_releff_FE
plot(pso_releff_FE, ref_line = 0)


# Predicted probabilities

pso_pred_FE <- <???>
pso_pred_FE
plot(pso_pred_FE, ref_line = c(0, 1))


## BONUS Q: Produce and plot the treatment ranks and rank probabilities for
##    each study in the network
##
##    (hint: use the posterior_ranks() and posterior_rank_probs() functions,
##    and remember that higher outcomes are better!)

# Ranks

pso_ranks_FE <- <???>
pso_ranks_FE
plot(pso_ranks_FE)


# Rank probabilities and cumulative rank probabilities

pso_rankprobs_FE <- <???>
pso_rankprobs_FE
plot(pso_rankprobs_FE)

pso_cumrankprobs_FE <- <???>
pso_cumrankprobs_FE
plot(pso_cumrankprobs_FE)



# Producing estimates for a specific target population --------------------

# Relative effects and ranks ----------------------------------------------

# Specify a data frame of mean covariate values for the target population
new_agd_means <- data.frame(
  bsa = 0.6,
  prevsys = 0.1,
  psa = 0.2,
  weight = 10,
  durnpso = 3)

## Q: Produce population-average treatment effects and rank probabilities for a
##    target population with the above mean covariate values

pso_releff_FE_new <- relative_effects(pso_fit_FE, <???>)
pso_releff_FE_new
plot(pso_releff_FE_new, ref_line = 0)

pso_rankprobs_FE_new <- posterior_rank_probs(pso_fit_FE, <???>,
                                             lower_better = FALSE)
pso_rankprobs_FE_new
plot(pso_rankprobs_FE_new)



# Predicted probabilities -------------------------------------------------

# For population-average event probabilities, we need to use numerical
# integration over the full joint covariate distribution.

new_agd_int <- data.frame(
  bsa_mean = 0.6,
  bsa_sd = 0.3,
  prevsys = 0.1,
  psa = 0.2,
  weight_mean = 10,
  weight_sd = 1,
  durnpso_mean = 3,
  durnpso_sd = 1
)

# Numerical integration for the joint distribution is set up in a similar manner
# to the network, except now we are adding the numerical integration points to
# the data frame `new_agd_int` which contains the covariate summaries.

new_agd_int <- add_integration(new_agd_int,
  durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
  prevsys = distr(qbern, prob = prevsys),
  bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
  weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
  psa = distr(qbern, prob = psa),
  cor = pso_net$int_cor,   # Use the weighted correlation matrix from the IPD studies
  n_int = 1000)

## Q: Complete the code below to produce population-average probabilities for a
##    population with the joint covariate distribution we have specified above.
##    Specify a N(-1.75, 0.08^2) distribution on the baseline log odds of response.

pso_pred_FE_new <- predict(pso_fit_FE, <???>)
pso_pred_FE_new
plot(pso_pred_FE_new, ref_line = c(0, 1))

