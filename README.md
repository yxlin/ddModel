# ddModel
ddModel provides density, distribution and random generation functions for the Diffusion Decision model (DDM). Supports model specification, parameter estimation, and likelihood computation.

# Getting Started

The package is mainly to support ggdmc, so you can use it together with other ggdmc supporting packages.

```
# Set up a stimuls drift rate model

model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
        t0 = "1", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0, s = 1, st0 = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "fastdm"
)

# Set up a population level distribution
pop_mean <- c(a = 1, sv = 0.1, sz = 0.25, t0 = 0.15, v = 2.5, z = 0.38)
pop_scale <- c(a = 0.05, sv = 0.01, sz = 0.01, t0 = 0.02, v = .5, z = 0.01)
pop_dist <- ggdmcPrior::BuildPrior(
    p0    = pop_mean,
    p1    = pop_scale,
    lower = c(0, 0, 0, 0, -10, 0),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

# Plot the population distribution 
ggdmcPrior::plot_prior(pop_dist)

# ---------------------------------------
# Subject-level and population-level response time models
sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

# Assume that a true parameter vector 
p_vector <- c(a = 1, sv = 0.1, sz = 0.25, t0 = 0.15, v = 2.5, z = .38)

# Use the simulate method in ddModel to simulate DDM model associated with the 
# `model` design:
# One subject
dat <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1)

# Multiple subjects
hdat <- simulate(pop_model, nsim = 128, n_subject = 32)
```

You couls also use the function in the package separately.

```
RT <- seq(0.1, 1.2, 0.01)
sz <- 0.05
sv <- 0.01
st0 <- 0.001

unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, sz = sz, sv = sv, t0 = 0.15,
    st0 = st0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
result <- pfastdm(RT, p_vector, is_lower = TRUE, debug = TRUE)

```


# Prerequisites
R (>= 3.3.0), ggdmcPrior, ggdmcModel, Rcpp (>= 1.0.7), RcppArmadillo (>= 0.10.7.5.0), and
ggdmcHeaders.

See DESCRIPTION for details

# Installation

From CRAN:
```
install.packages("ddModel")
```