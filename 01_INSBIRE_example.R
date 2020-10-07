#### INSBIRE ####
## Interactions in Nutrient Spirals using BayesIan REgression
## Author: Matthias Pucher
## matthias.pucher@boku.ac.at
## This example script was developed during the work on
## Pucher et al. (2020): Complex interactions of in-stream
## DOM and nutrient spiralling unravelled by Bayesian
## regression analysis.

#### Packages: install and load ####
## Please see the vignettes and documentation of these packages for
## further information.

#install.packages("tidyverse","brms","shinystan")

require(tidyverse)
require(brms)

#### Data: generate artificial data ####
## We generate data following the nutrient spiralling concept.
## The data generation is artificial, but should be able to demonstrate
## the INSBIRE approach in R.
## The user can change our assumptions or add noise to see the influences
## on the final results.

## for reproducible random data, we set a seed
## for a different random realisation 
set.seed(42)

## C3, Concentration of substance 3, is generated randomly for
## 10 experiment runs and 12 sequential sampling points

C3 <- rnorm(10, 10, 2) %>% # 10 runs, mean of 10, standard deviation of 2
  sapply(function(m){ rnorm(12, m, 0.1) }) # 12 samplings per run, mean from distribution above, sd of 0.1

## C2 is correlated to C3, just additional noise is added.
C2 <- C3 + rnorm(length(C3),0,1.3) # add noise with a mean of 0 and a sd of 1.3

## check the correlation
cor(as.vector(C2),as.vector(C3))

## C4 is generated indepently from the other substances
C4 <- rnorm(ncol(C3), 3, 1) %>% # same dimension as C3, mean of 3, sd of 1
  sapply(function(m){ rnorm(nrow(C3), m, 0.1) }) # slight deviation (sd of 1) from the mean of the run

## assuming values for the parameters
k <- 4.7
m3 <- 0.5
m4 <- -0.3

## vf is constructed, optionally noise can (and should) be added
vf <- (k * C3 ^ m3 * C4 ^ m4) + rnorm(length(C3), 0, 0.05)

## C1, concentration of solvent 1, this will be a dependent variable
## we generate only the data of the first sampling, the leachate introduction
C1 <- matrix(NA, nrow = nrow(C3), ncol = ncol(C3))
C1[1,] <- runif(ncol(C3), 2000, 3500) # uniformly distributed values between 2000 and 3500 as starting conditions

## random ambient concentration
C1_amb <- matrix(rnorm(length(C3), 1800, 30), nrow = nrow(C3), ncol = ncol(C3)) # ambient concentrations are around 1800

## sampling points have a random distance to each other, similar for all runs
stretch <- matrix(c(NA, runif(nrow(C3) - 1, 25, 30)), nrow = nrow(C3), ncol = ncol(C3)) # the stretches have a length between 25 and 30 m

## u*z, flow velocity * water depth, is generated for each run
uz <- matrix(runif(ncol(C3), 0.03, 0.06), nrow = nrow(C3), ncol = ncol(C3), byrow = TRUE) # small stream with u*z uniformly distributed between 0.03 and 0.06

## sequential calculation of the retention of C1
for(col in 1:ncol(C1)){
  for(row in 2:nrow(C1)){
    C1[row, col] <- C1_amb[row, col] + (C1[row - 1, col] - C1_amb[row, col]) * exp(- stretch[row, col] * vf [row, col] / (uz[row, col] * 1000 * 60 )) # 1000 * 60 added to get mm min¯¹ for the uptake velocity
  }
}

C1

C1 <- C1 + rnorm(length(C1), 0, 5) # noise can be added optionally

## combine data to one data.frame
data <- lapply(list(C1 = C1, C2 = C2, C3 = C3, C4 = C4, C1_amb = C1_amb, stretch = stretch, uz = uz), as.vector) %>%
  bind_cols() %>%
  mutate(C1_lag = lag(C1)) %>%
  mutate(D = C1_lag - C1)

#### INSBIRE ####
## In this section, we allpy the INSBIRE approach on the artificially
## created data. The user can use his or her own data from here.
## BEWARE: the calculation takes several minutes to hours depending on the computer, the data, the formula and the selected priors!

#### vf ####
## vf is calculated as a parameter
## the formula is Eq. (4) from the paper
model_vf <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * vf / (uz * 1000 * 60))),
                             vf ~ 1,
                             nl = TRUE),
                data = data,
                family = gaussian(),
                prior = c(
                  prior(lognormal(log(2), 3), nlpar = "vf", lb = 0.01, ub = 35)),
                cores = 4,
                save_all_pars = TRUE,
                sample_prior = "yes"
)

summary(model_vf)

#### interactions ####
## Equation (8) of the paper
## include interactions with C2
model_C2 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C2^m2 / (uz * 1000 * 60))),
                             k ~ 1,
                             m2 ~ 1,
                             nl = TRUE),
                data = data,
                family = gaussian(),
                prior = c(
                  prior(normal(0.2,0.4), nlpar = "m2", lb = -1, ub = 1),
                  prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                cores = 4,
                save_all_pars = TRUE,
                sample_prior = "yes"
)

pairs(model_C2)
summary(model_C2)

## how does the model including C2 perform incomparison to no additional information for vf
bf_vf_C2 <- bayes_factor(model_vf, model_C2)
bf_vf_C2

# the model with C2-interactions included performs much better
# we ask for the favour of the vf-model, so the inverse is needed here.
1/bf_vf_C2$bf
# it is around 1000 times more likely. This would mean a decisive improvement over the simpler model

model_C3 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C3^m3 / (uz * 1000 * 60))),
                             k ~ 1,
                             m3 ~ 1,
                             nl = TRUE),
                data = data,
                family = gaussian(),
                prior = c(#prior(gamma(2, .1), class = nu),
                  prior(normal(0.2,0.4), nlpar = "m3", lb = -1, ub = 1),
                  prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                control = list(adapt_delta = 0.95, max_treedepth = 20),
                cores = 4,
                save_all_pars = TRUE,
                sample_prior = "yes"
)

pairs(model_C3)
summary(model_C3)

bf_vf_C3 <- bayes_factor(model_vf, model_C3)
bf_vf_C3

# the model with C3-interactions included performs much better
# we ask for the favour of the vf-model, so the inverse is needed here.
1/bf_vf_C3$bf
# it is around 10^10 times more likely. This would mean a decisive improvement over the simpler model

# we can directly see a higher Bayes factor than with model_C2, but we can also compare the two with each other
bf_C2_C3 <- bayes_factor(model_C2, model_C3)
bf_C2_C3
1/bf_C2_C3$bf

# model_C3 is around 10^7 times more likely, so again an improvement
# since we added more noise to C2 than to C3 during the data construction, this is reasonable

model_C4 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C4^m4 / (uz * 1000 * 60))),
                             k ~ 1,
                             m4 ~ 1,
                             nl = TRUE),
                data = data,
                family = gaussian(),
                prior = c(#prior(gamma(2, .1), class = nu),
                  prior(normal(-0.2,0.4), nlpar = "m4", lb = -1, ub = 1),
                  prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                #control = control[[response]][[mod]],
                #iter = iter[[response]][[mod]],
                cores = 4,
                #chains = chains[[response]][[mod]],
                save_all_pars = TRUE,
                sample_prior = "yes"
)

pairs(model_C4)
summary(model_C4)

bf_vf_C4 <- bayes_factor(model_vf, model_C4)
bf_vf_C4

# the model with C4-interactions included performs much better
# we ask for the favour of the vf-model, so the inverse is needed here.
1/bf_vf_C4$bf
# it is around 10^16 times more likely. This would mean a decisive improvement over the simpler model

# now we combine different interactions
# adapt_delta is increased, because the default value of 0.8 was not enough for the model to converge
model_C2C3 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C2^m2 * C3^m3 / (uz * 1000 * 60))),
                             k ~ 1,
                             m2 ~ 1,
                             m3 ~ 1,
                             nl = TRUE),
                data = data,
                family = gaussian(),
                prior = c(#prior(gamma(2, .1), class = nu),
                  prior(normal(0.2,0.4), nlpar = "m2", lb = -1, ub = 1),
                  prior(normal(0.2,0.4), nlpar = "m3", lb = -1, ub = 1),
                  prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                control = list(adapt_delta = 0.90),
                cores = 4,
                save_all_pars = TRUE,
                sample_prior = "yes"
)

pairs(model_C2C3)
summary(model_C2C3)

# can the model be improved if C2 and C3 are used together?
bf_C3_C2C3 <- bayes_factor(model_C3,model_C2C3)
bf_C3_C2C3

# there is slight favour for model_C3, so we would reject model_C2C3. Even with a slight favour for model_C2C3 it can be argued to reject it due to its higher complexity.

model_C3C4 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C3^m3 * C4^m4 / (uz * 1000 * 60))),
                               k ~ 1,
                               m3 ~ 1,
                               m4 ~ 1,
                               nl = TRUE),
                  data = data,
                  family = gaussian(),
                  prior = c(#prior(gamma(2, .1), class = nu),
                    prior(normal(0.2,0.4), nlpar = "m3", lb = -1, ub = 1),
                    prior(normal(-0.2,0.4), nlpar = "m4", lb = -1, ub = 1),
                    prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                  iter = 2000,
                  cores = 4,
                  save_all_pars = TRUE,
                  sample_prior = "yes"
)

pairs(model_C3C4)
summary(model_C3C4)

# can the model be improved by adding both C3 and C4?
bf_C3_C3C4 <- bayes_factor(model_C3,model_C3C4)
bf_C3_C3C4
1/bf_C3_C3C4$bf

# Yes, the combnation of both components is 10^12 times more likely.
# the coefficients are approximately k = 5.4, m3 = 0.45 and m4 = -0.31
# a comparison with the original values (k = 4.7, m3 = 0.5, m4 = -0.3) shows,
# that we could determine the basic behaviour of the relations

model_C2C3C4 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C2 ^ m2 * C3^m3 * C4^m4 / (uz * 1000 * 60))),
                               k ~ 1,
                               m2 ~ 1,
                               m3 ~ 1,
                               m4 ~ 1,
                               nl = TRUE),
                  data = data,
                  family = gaussian(),
                  prior = c(#prior(gamma(2, .1), class = nu),
                    prior(normal(0.2,0.4), nlpar = "m2", lb = -1, ub = 1),
                    prior(normal(0.2,0.4), nlpar = "m3", lb = -1, ub = 1),
                    prior(normal(-0.2,0.4), nlpar = "m4", lb = -1, ub = 1),
                    prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                  control = list(adapt_delta = 0.95, max_treedepth = 20),
                  cores = 4,
                  save_all_pars = TRUE,
                  sample_prior = "yes"
)

pairs(model_C2C3C4)
summary(model_C2C3C4)

bf_C3C4_C2C3C4 <- bayes_factor(model_C3C4,model_C2C3C4)
bf_C3C4_C2C3C4

# the simpler model,including C3 and C4 only is 2 times more likely. The coefficients are still close to the original values. m2 is most probable around 0, which seems reasonable since the factor has no influence on the model then.
# model_C3C4 would be our best guess.

#### recalculate model with more iterations ####
model_C3C4 <- brm(formula = bf(D ~ -1*(-C1_lag + C1_amb + (C1_lag - C1_amb) * exp(- stretch * k * C3^m3 * C4^m4 / (uz * 1000 * 60))),
                               k ~ 1,
                               m3 ~ 1,
                               m4 ~ 1,
                               nl = TRUE),
                  data = data,
                  family = gaussian(),
                  prior = c(#prior(gamma(2, .1), class = nu),
                    prior(normal(0.2,0.4), nlpar = "m3", lb = -1, ub = 1),
                    prior(normal(-0.2,0.4), nlpar = "m4", lb = -1, ub = 1),
                    prior(lognormal(log(2), 3), nlpar = "k", lb = 0.01, ub = 40)),
                  control = list(adapt_delta = 0.95, max_treedepth = 20),
                  iter = 10000,
                  cores = 8,
                  chains = 8,
                  save_all_pars = TRUE,
                  sample_prior = "yes"
)

#### plot posterior ####
## diagnose and explore the model
shinystan::launch_shinystan(model_C3C4)

#### plot relations ####
simData <- model_C3C4 %>%
  posterior_samples()

## extract and restructure data for the plot
simData2 <- lapply(c("C3","C4") %>% setNames(.,.), function(par){ # create a data grid for C3 and C4 within the observed range
  range <- range(data[[par]])
  seq <- seq(from = min(data[[par]], na.rm = TRUE), to = max(data[[par]], na.rm = TRUE), length.out = 200)
}) %>%
  expand.grid() %>%
  slice_sample(n  = nrow(simData), replace = TRUE) %>% # repeat C3 and C4 to the number of simulated posterior samples
  bind_cols(simData, .) %>%
  mutate(response = "C1") %>%
  mutate(vf = b_k_Intercept * C3 ^ b_m3_Intercept * C4 ^ b_m4_Intercept) %>% # calculate vf from the simulated posterior samples and the concentrations
  pivot_longer(one_of("C3", "C4"), names_to = "parameter", values_to = "conc") %>%
  group_by(parameter, conc) %>%
  mutate(perc = ifelse(vf > quantile(vf,0.495) & vf < quantile(vf, 0.505), "50", ifelse(vf > quantile(vf, 0.25) & vf < quantile(vf, 0.75), "50", ifelse(vf > quantile(vf, 0.05) & vf < quantile(vf, 0.95), "90", NA)))) %>% # distinguish between quantiles for the according colour in the plot
  ungroup() %>%
  arrange(conc)

# plot the data
simData2 %>%
  filter(!is.na(perc)) %>% # filter points beyond the 5 and 95 percentiles
  ggplot(aes(x = conc, y = vf))+
  geom_jitter(aes(colour = perc), size = 1, alpha = 0.2)+
  facet_wrap(parameter~., scales= "free", strip.position = "bottom")+
  labs(y ="uptake length vf in mm min¯¹" ,x="")
