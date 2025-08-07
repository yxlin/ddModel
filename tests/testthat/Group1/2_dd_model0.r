#  q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ddModel", "ggdmcPrior", "ggdmcPhi")
suppressPackageStartupMessages(tmp <- sapply(pkg, require,
    character.only = TRUE
))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ddModel/tests/testthat/Group1/"
model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
        t0 = "1", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0.25, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("r1", "r2"), # r1 = lower; r2 = upper
    type = "fastdm"
)

pop_mean <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = 0.38)
pop_scale <- c(a = 0.05, sz = 0.01, t0 = 0.02, v = .5, z = 0.01)
pop_dist <- ggdmcPrior::BuildPrior(
    p0    = pop_mean,
    p1    = pop_scale,
    lower = c(0, 0, 0, -5, 0),
    upper = rep(10, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)


sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

unsorted_p_vector <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = .38)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]

# p_vector <- c(A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5, t0 = .15)

# trials <- simulate_ddm_trials(
#     rt_model_r = sub_model,
#     parameters_r = p_vector,
#     n_trial = as.integer(20),
#     debug = TRUE
# )
# 0 = lower; 1 = upper


dat <- simulate(sub_model, nsim = 32, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 64, n_subject = 5)

# result <- .new_convert2datalist(dat)
# result[[1]]
# result[[2]]

# result <- .new_convert2datalist(hdat)
# result[[1]]
# result[[2]]

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
sub_dmis2 <- ggdmcModel::BuildDMI(hdat, model)

sub_dmis[[1]]@data
sub_dmis[[1]]@node_1_index


sub_dmis2[[1]]@data
sub_dmis2[[1]]@node_1_index

sub_dmis2[[2]]@data
sub_dmis2[[2]]@node_1_index

sub_dmis2[[3]]@data
sub_dmis2[[3]]@node_1_index

order(names(sub_dmis2[[1]]@data))
sub_dmis2[[2]]
head(trials)
trials
