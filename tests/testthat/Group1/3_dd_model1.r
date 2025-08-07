q(save = "no")

cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ddModel", "ggdmcPrior", "ggdmcPhi")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ddModel/tests/testthat/Group1/"

model <- ggdmcModel::BuildModel(
    p_map = list(a = "1", v = "S", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1", s = "1", precision = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0.25, s = 1, st0 = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "fastdm"
)

hyper_model <- ggdmcModel::BuildModel(
    p_map = list(a = "1", v = "S", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1", s = "1", precision = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0.25, s = 1, st0 = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "hyper"
)


pop_mean <- c(a = 2, sv = 1, sz = 0.3, t0 = 0.3, v.s1 = 4, v.s2 = 3, z = 0.5)
pop_scale <- c(a = 0.5, sv = 0.3, sz = 0.1, t0 = 0.05, v.s1 = .5, v.s2 = .5, z = 0.1)
pop_dist <- ggdmcPrior::BuildPrior(
    p0    = pop_mean,
    p1    = pop_scale,
    lower = rep(0, 0, 0, 0, -5, -5, 0),
    upper = rep(10, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

# ggdmcPrior::plot_prior(pop_dist)

sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

p_vector <- c(a = 1, sv = 0.2, sz = 0.25, t0 = 0.15, v.s1 = 2.2, v.s2 = 2, z = .38)

trials <- simulate_ddm_trials(
    rt_model_r = sub_model,
    parameters_r = p_vector,
    n_trial = as.integer(20),
    debug = TRUE
)
# 0 = lower; 1 = upper


dat <- simulate(sub_model, nsim = 32, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 32, n_subject = 3)

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)

sub_dmis[[1]]
mean(dat$C)
mean(hdat$C)
ps <- attr(hdat, "parameters")
ps

true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)


# Generate subj samples -------------------------
p0 <- rep(0, model@npar)
names(p0) <- model@pnames

p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)

nmc <- 2
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames)
print(sub_theta_input)



# sum(sapply(likelihoods, ggdmcLikelihood:::.sumlog))

# dmi <- sub_dmis[[1]]
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 846671, verbose = TRUE)
sub_samples
sub_dmis[[1]]

# save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
#     sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
#     file = save_path
# )


# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = c(0, 0, 0, 0, 0),
    upper = rep(NA, model@npar),
    dist = rep("tnorm", model@npar),
    log_p = rep(TRUE, model@npar)
)

# Prior log likelihoods
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
l_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
s_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
nmc
pop_priors <- ggdmcPrior::set_priors(p_prior = model_likelihood, l_prior = l_prior, s_prior = s_prior)
pop_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = pop_priors@pnames)

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)
pop_samples[[1]]
