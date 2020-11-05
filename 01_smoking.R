## Example: Smoking Cessation

library(multinma)
options(mc.cores = parallel::detectCores())



# Data --------------------------------------------------------------------

head(smoking)


# Network -----------------------------------------------------------------

smknet <- set_agd_arm(smoking,
                      study = studyn,
                      trt = trtc,
                      r = r,
                      n = n,
                      trt_ref = "No intervention")
smknet

plot(smknet, weight_edges = TRUE, weight_nodes = TRUE)


# FE NMA ------------------------------------------------------------------

summary(normal(scale = 100))

smkfit_FE <- nma(smknet,
                 trt_effects = "fixed",
                 prior_intercept = normal(scale = 100),
                 prior_trt = normal(scale = 100))

smkfit_FE


# Plot prior vs posterior

plot_prior_posterior(smkfit_FE)


# Check DIC and residual deviance

smk_dic_FE <- dic(smkfit_FE)
smk_dic_FE

plot(smk_dic_FE)


# RE NMA ------------------------------------------------------------------

summary(normal(scale = 100))
summary(half_normal(scale = 5))

## Q: Modify the code to fit a random effects model, with a half-N(0, 5^2)
##    prior on tau

## Tip: To get help on a function use `?`
##      E.g. run ?nma for help on the nma() function

smkfit_RE <- nma(smknet,
                 trt_effects = "fixed",
                 prior_intercept = normal(scale = 100),
                 prior_trt = normal(scale = 100),
                 prior_het = <???>)

smkfit_RE


# Plot prior vs posterior for tau
plot_prior_posterior(smkfit_RE, prior = "het")

## Q: Check the DIC and plot the residual deviance contributions.
##    Compare with the FE model.

# Calculate DIC
smk_dic_RE <- <???>

# Compare to FE model
smk_dic_RE
smk_dic_FE

# Check residual deviance contributions
plot(smk_dic_RE)

smoking[smoking$r == 0, ]


# Fit UME model -----------------------------------------------------------

## Q: Modify the code to fit a RE UME model to check for inconsistency.
##    Compare DIC with the consistency RE model, and produce a dev-dev plot.

smkfit_UME <- nma(smknet,
                  consistency = "consistency",
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = half_normal(scale = 5))
smkfit_UME


# Compare DIC and residual deviance
smk_dic_RE

smk_dic_UME <- dic(smkfit_UME)
smk_dic_UME

# Produce dev-dev plot
plot(smk_dic_RE, smk_dic_UME, show_uncertainty = FALSE)
plot(smk_dic_RE, smk_dic_UME, interval_alpha = 0.1)


# Relative effects --------------------------------------------------------

## Q: Produce relative effect estimates for all treatment contrasts, from the
##    RE model.

smk_releff <- relative_effects(<???>)
smk_releff
plot(smk_releff, ref_line = 0)

## BONUS Q: Transform these relative effects into odds ratios
##          (hint: use as.array() and summary())



# Absolute effects --------------------------------------------------------

## Q: Produce predicted probabilities of cessation for a population with
##    baseline log odds (for no intervention) distributed as N(-2.34, 0.12^2)

smk_pred <- predict(smkfit_RE,
                    baseline = <???>,
                    trt_ref = "No intervention",
                    type = <???>)

smk_pred
plot(smk_pred, ref_line = c(0, 1))


## BONUS Q: Investigate the different plotting options. Can you produce a density
##          plot of the predicted probabilities?



# Ranks -------------------------------------------------------------------

## Q: Produce posterior ranks of each treatment (hint: higher treatment effects
##    are better - greater odds of cessation)

smk_ranks <- posterior_ranks(<???>)
smk_ranks
plot(smk_ranks)


# Rank probabilities ------------------------------------------------------

## Q: Produce posterior rank probabilities and cumulative rank probabilities

# Rank probabilities
smk_rankprobs <- posterior_rank_probs(<???>)
smk_rankprobs
plot(smk_rankprobs)

# Cumulative rank probabilities
smk_cumrankprobs <- posterior_rank_probs(<???>, cumulative = TRUE)
smk_cumrankprobs
plot(smk_cumrankprobs)

# We can customise the plots using ggplot commands
plot(smk_cumrankprobs) +
  ggplot2::facet_null() +
  ggplot2::aes(colour = Treatment, linetype = Treatment)
