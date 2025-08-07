# 📦 ddModel: Core Functions for the Decision Diffusion Model

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ddModel)](https://cran.r-project.org/package=ddModel)
[![Downloads](https://cranlogs.r-pkg.org/badges/ddModel)](https://cran.r-project.org/package=ddModel)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/yxlin/ddModel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yxlin/ddModel/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://codecov.io/gh/yxlin/ddModel/branch/main/graph/badge.svg)](https://codecov.io/gh/yxlin/ddModel)

`ddModel` provides fast and flexible computational tools for the **Decision Diffusion Model (DDM)**, a cognitive model used to jointly analyse *choice* and *response time* (RT) data in speeded decision-making tasks.

## 🔍 Overview
This package supports:

- 📈 Density, distribution, and random sampling functions for the DDM
- ⚙️ Flexible parameter specification, allowing you to:
  - Fix parameters globally
  - Constrain them by experimental conditions
  - Vary them trial-by-trial

- 🧪 Simulation tools for generating synthetic datasets under known parameters
- 📊 Likelihood evaluation functions for integration with external optimisation or MCMC routines

>❗ Note: ddModel, although does not include an optimiser, is designed to integrate with the `ggdmc` package, enabling the Bayesian parameter optimisation. Nevertheless, the user may use other estimation pipelines.

## 🧠 Features
- Implements canonical DDM components:
`a` (boundary separation), `v` (drift rate), `t₀` (non-decision time), `z` (starting point), plus optional variability parameters `sv`, `sz`, `st₀`
- Fully vectorised functions for efficient simulation and likelihood evaluation on large datasets
- No dependencies beyond base R and Rcpp; integrates easily with the ggdmc ecosystem

## 📌 Typical Use Cases
`ddModel` is useful for:
- Experimental psychologists modelling two-alternative forced-choice (2AFC) tasks
- Cognitive scientists conducting simulation-based model recovery or parameter estimation
- Researchers generating synthetic data for power analyses or benchmarking

## 🚀 Getting Started
While `ddModel` can be used independently, it is designed to work with the `ggdmc` ecosystem, particularly:
- ggdmcModel
- ggdmcPrior
- ggdmcHeaders

## ✅ Example Workflow


```r
# Load packages
library(ggdmcModel)
library(ggdmcPrior)
library(ddModel)

# Set up a stimulus drift rate model
model <- BuildModel(
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

# Set up a population-level prior distribution
pop_mean  <- c(a = 1, sv = 0.1, sz = 0.25, t0 = 0.15, v = 2.5, z = 0.38)
pop_scale <- c(a = 0.05, sv = 0.01, sz = 0.01, t0 = 0.02, v = 0.5, z = 0.01)
pop_dist  <- BuildPrior(
  p0    = pop_mean,
  p1    = pop_scale,
  lower = c(0, 0, 0, 0, -10, 0),
  upper = rep(NA, length(pop_mean)),
  dists = rep("tnorm", length(pop_mean)),
  log_p = rep(FALSE, length(pop_mean))
)

# Visualise the prior
plot_prior(pop_dist)

# Subject-level and population-level model setup
sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

# Simulate subject-level data
p_vector <- c(a = 1, sv = 0.1, sz = 0.25, t0 = 0.15, v = 2.5, z = 0.38)
dat      <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1)

# Simulate hierarchical data (32 subjects)
hdat     <- simulate(pop_model, nsim = 128, n_subject = 32)
```

## ⚙️ Core Function Example (using `pfastdm`)


```r

RT <- seq(0.1, 1.2, 0.01)
params <- c(
  a = 1, v = 1.5, zr = 0.5, d = 0,
  sz = 0.05, sv = 0.01, t0 = 0.15, st0 = 0.001,
  s = 1, precision = 3
)
# Ensure parameter names are ordered
params <- params[sort(names(params))]

# Compute lower-bound response density
result <- pfastdm(RT, params, is_lower = TRUE, debug = TRUE)

```

## 🧩 Dependencies
- R (≥ 3.3.0)
- `Rcpp` (≥ 1.0.7)
- `RcppArmadillo` (≥ 0.10.7.5.0)
- `ggdmcModel`, `ggdmcPrior`, `ggdmcHeaders`

## 📦 Installation
Install the latest release from CRAN:


```r
install.packages("ddModel")
```

Or install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("yxlin/ddModel")
```


## 📚 References
- Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the parameters of the diffusion model: An empirical validation. *Behavior Research Methods, Instruments, & Computers*, 36(3), 347–360. [https://doi.org/10.3758/BF03196893](https://doi.org/10.3758/BF03196893)
- Voss, A., & Voss, J. (2007). Fast-dm: A free program for efficient diffusion model analysis. *Behavior Research Methods*, 39(4), 767–775. [https://doi.org/10.3758/BF03192967](https://doi.org/10.3758/BF03192967)
- Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for two-choice decision tasks. *Neural Computation*, 20(4), 873–922. [https://doi.org/10.1162/neco.2008.12-06-420](https://doi.org/10.1162/neco.2008.12-06-420)


## 📬 Contact
For questions, feedback, or contributions, please contact:

- **Yi-Shin Lin**
- Email: [yishinlin001@gmail.com](mailto:yishinlin001@gmail.com)
- GitHub: @yxlin


