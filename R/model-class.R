#' @importClassesFrom ggdmcModel model
#' @importClassesFrom ggdmcPrior prior
NULL

#' An S4 Class Representing a Decision Diffusion Model
#'
#' @description
#' The \code{ddm} class represents a complete Decision Diffusion Model
#' (DDM) specification. It encapsulates a \code{model} objects from
#' the \pkg{ggdmcModel} package, and optionally a population distribution
#' created via the \code{BuildPrior} function from the \pkg{ggdmcPrior}
#' package.
#'
#' @slot model A \code{model} object from \pkg{ggdmcModel} package
#' that defines the DDM structure, including parameter names, mapping,
#' and data attributes.
#' @slot population_distribution An optional prior distribution, typically
#' created via \code{BuildPrior} from the \pkg{ggdmcPrior} package. This
#' slot is primarily used for parameters or hierarchical modelling and can
#' be \code{NULL}.
#'
#' @param model A required \code{model} object defining the model
#' specification.
#' @param population_distribution Optional A population-level prior
#' distribution, typically constructed with \code{BuildPrior}. This argument
#' enables simulations of subject-level parameters. Default is \code{NULL}.
#'
#' @return An object of S4 class \code{ddm}, which includes:
#' \itemize{
#' \item \code{model}: the model specification
#' \item \code{population_distribution}: the prior distribution (if provided)
#' }
#'
#' @details
#' The constructor function \code{setDDM} creates a \code{ddm} object by
#' wrapping the model and optionally a population distribution into a single
#' structure. This class is designed for compatibility with sampling and
#' simulation functions within the \pkg{ggdmc} ecosystem.
#'
#' @seealso \code{\link[ggdmcModel]{model-class}}, \code{\link[ggdmcPrior]{BuildPrior}}
#'
#' @rdname ddm-class
#' @export
setClass("ddm",
  slots = list(
    model = "model",
    population_distribution = "ANY"
  ),
  validity = function(object) {
    if (!methods::is(object@model, "model")) {
      return("Slot 'model' must be a 'model' object from ggdmcModel")
    }
    if (!is.null(object@population_distribution) &&
      !methods::is(object@population_distribution, "list")) {
      return("population_distribution must be NULL or a list.")
    }
    TRUE
  }
)

#' @importFrom methods new
#' @rdname ddm-class
#' @export
setDDM <- function(model, population_distribution = NULL) {
  out <- new("ddm",
    model = model,
    population_distribution = population_distribution
  )
  out
}

#' @importFrom ggdmcPrior rprior
.validate_parameters <- function(
    parameter_matrix, pop_dist, rt_model,
    max_attempts = 1000) {
  # Validate each subject's parameters, resampling invalid ones
  # Returns validated matrix or NULL if max attempts reached
  n_subject <- nrow(parameter_matrix)

  all_valid <- rep(FALSE, n_subject)
  for (i in seq_len(n_subject)) {
    current_params <- parameter_matrix[i, ]
    all_valid[i] <- validate_ddm_parameters(rt_model, current_params)
  }

  for (i in seq_len(n_subject)) {
    if (!all_valid[i]) {
      message("Subject ", i, " has invalid parameters. Resample a new set from the population distribution")

      for (attempt in seq_len(max_attempts)) {
        tmp_params <- rprior(pop_dist, n = 1)
        is_valid <- validate_ddm_parameters(rt_model, tmp_params[, 1])
        if (is_valid) {
          cat("After", attempt, "attempts, I found a new set of parameters for Subject", i, "\n")

          parameter_matrix[i, ] <- tmp_params[, 1]
          all_valid[i] <- is_valid
          break
        }
      }
    }
  }

  if (!all(all_valid)) {
    stop("Failed to generate valid parameters for all subjects after ", max_attempts, " attempts")
  }

  return(parameter_matrix)
}

#' @importFrom ggdmcPrior rprior
.prepare_parameter_matrix <- function(
    model, n_subject, pop_dist = NULL,
    seed = NULL) {
  # same function as in lbaModel

  parameter_names <- model@pnames

  # Validate input arguments
  if (is.null(pop_dist)) {
    stop("You must provide a population distribution.")
  }

  # Validate all model parameters have corresponding priors
  if (!all(parameter_names %in% names(pop_dist))) {
    stop("Model parameters differ from the population distribution.")
  }

  # Sample parameters for each subject from the priors
  if (is.null(seed)) {
    seed <- as.integer(Sys.time()) %% 1000000 # Generate a random seed
    message("No seed provided. Using R-generated seed: ", seed)
    set.seed(seed)
  } else {
    set.seed(seed)
  }

  t(ggdmcPrior::rprior(pop_dist, n = n_subject))
}

.separate_condition_column <- function(
    data, factors, responses = NULL,
    match_map = NULL) {
  # Search "DDM" for differenc
  # Split the "Condition" column into parts based on the delimiter "."
  condition_parts <- strsplit(data$Condition, "\\.")

  # Create a new data frame to store the separated columns
  separated_data <- data.frame(data)

  # Iterate over the factors and extract the corresponding parts
  for (factor_name in names(factors)) {
    # Extract the part of the condition corresponding to the current factor
    separated_data[[factor_name]] <- sapply(condition_parts, function(x) {
      # Find the index of the factor in the condition parts
      index <- which(x %in% factors[[factor_name]])
      if (length(index) == 0) {
        return(NA) # Return NA if the factor is not found
      }
      return(x[index])
    })
  }

  # Drop the original "Condition" column
  separated_data <- separated_data[, !names(separated_data) %in% "Condition"]

  # DDM assumes hitting lower = FALSE, and upper = TRUE
  separated_data$C <- ifelse(separated_data$R == 0, FALSE, TRUE)
  if (!is.null(responses) && !is.null(match_map)) {
    separated_data$R <- ifelse(
      separated_data$C,
      sapply(separated_data$S, function(s) match_map$M[[s]]),
      sapply(separated_data$S, function(s) setdiff(responses, match_map$M[[s]]))
    )
  }

  # DDM convert additonal C columns.
  columns_to_skip <- c("RT", "C") # Define columns to exclude
  for (col in names(separated_data)) {
    if (!col %in% columns_to_skip) {
      separated_data[[col]] <- as.factor(separated_data[[col]])
    }
  }

  return(separated_data)
}

#' Simulate Design-based Decision Diffusion Model (DDM) Trials
#'
#' Generates synthetic response time (RT) and choice data using a diffusion
#' decision model (DDM) with specified parameters. The simulation allows for
#' flexible parameter settings and time-domain controls.
#'
#' @name simulate_ddm_trials
#'
#' @param rt_model_r An S4 object of class \code{ddm} specifying the DDM
#' configuration. Must contain slots for boundary separation, drift rate,
#' and other DDM-relevant parameters.
#'
#' @param parameters_r A numeric vector of DDM parameters. Expected elements:
#'    \describe{
#'      \item{a}{Boundary separation (threshold)}
#'      \item{v}{Drift rate (evidence accumulation speed)}
#'      \item{t0}{Non-decision time (encoding + motor response baseline)}
#'      \item{z}{Starting point (bias; must satisfy \code{0 < z < a})}
#'      \item{d}{(Optional) Difference in motor execution time between
#' boundaries}
#'      \item{sv}{(Optional) Inter-trial variability in drift rate}
#'      \item{st0}{(Optional) Inter-trial variability in non-decision time}
#'      \item{sz}{(Optional) Inter-trial variability in starting point}
#'    }
#'
#' @param time_parameters_r A numeric vector controlling simulation's
#' temporal dynamics:
#'        \describe{
#'          \item{[1]}{Minimum time (default: \code{-0.5})}
#'          \item{[2]}{Maximum time (default: \code{0.5})}
#'          \item{[3]}{Time step for numerical integration (default:
#' \code{0.01})}
#'        }
#'
#' @param n_trial Integer. Number of trials to simulate (default: \code{1}).
#'
#' @param debug Logical. If \code{TRUE}, prints internal simulation states
#' and timing diagnostics.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{rt}{Simulated response times in seconds}
#'   \item{response}{Binary decision outcome
#' (0 = lower boundary, 1 = upper boundary)}
#' }
#'
#' @details
#' The core simulation is implemented in C++ via \code{simulate_ddm_trials}
#' for high performance in large-scale experiments. A complementary R
#' wrapper (\code{.simulate_ddm_trials}) provides additional functionality,
#' such as random seed control and input validation.
#'
#' The function \code{validate_ddm_parameters} checks whether the input
#' parameter vector satisfies the constraints of the diffusion model.
#'
#' @references
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory
#' and data for two-choice decision tasks. Neural Computation, 20(4), 873-922.
#'
#' @seealso \code{\link{ddm-class}}, \code{\link{setDDM}}, \code{\link{validate_ddm_parameters}}
#'
#' @examples
#' # Basic simulation with default parameters
#' if (requireNamespace("ggdmcModel", quietly = TRUE)) {
#'   BuildModel <- getFromNamespace("BuildModel", "ggdmcModel")
#'   model <- ggdmcModel::BuildModel(
#'     p_map = list(
#'       a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
#'       t0 = "1", st0 = "1", s = "1", precision = "1"
#'     ),
#'     match_map = list(M = list(s1 = "r1", s2 = "r2")),
#'     factors = list(S = c("s1", "s2")),
#'     constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
#'     accumulators = c("r1", "r2"),
#'     type = "fastdm"
#'   )
#' }
#' sub_model <- setDDM(model)
#' p_vector <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = .38)
#' time_parameters <- c(t_min = -0.5, tmax = 0.5, dt = 0.01)
#'
#' # Simulation with custom time parameters
#' sim_data <- simulate_ddm_trials(
#'   rt_model_r = sub_model,
#'   parameters_r = p_vector,
#'   time_parameters_r = time_parameters,
#'   n_trial = 32
#' )
#'
#' @export
NULL


.simulate_ddm_trials <- function(
    rt_model_r, parameters_r, time_parameters_r,
    n_trial, seed, debug) {
  if (debug) {
    message("Set seed to: ", seed)
  }
  set.seed(seed)
  trials <- simulate_ddm_trials(
    rt_model_r = rt_model_r,
    parameters_r = parameters_r,
    time_parameters_r = time_parameters_r,
    n_trial = n_trial,
    debug = debug
  )
  return(trials)
}



#' Simulate Data from a Decision Diffusion Model
#'
#' Simulate response times and choices from a Decision Diffusion model (DDM)
#' associated with a specific experimental design. The specification is
#' typically specified, using \code{ggdmcModel::BuildModel}).
#'
#' @param object An object of class \code{"ddm"} that defines the model
#' structure and parameters.
#' @param nsim Integer. Number of trials to simulate per subject. Defaults
#' to \code{4}. The function adjusts \code{nsim} if not divisible by
#' number of conditions.
#' @param seed Optional integer. Sets the random seed for reproducibility.
#' Defaults to \code{NULL}.
#' @param n_subject Integer. Number of subjects to simulate. Defaults to
#' \code{3}.
#' @param parameter_vector Optional. A named vector or list of parameters
#' (e.g., \code{a}, \code{z}, \code{v}, \code{t0}) used as the true value
#' to simulate. The function by default generates one or many true
#' parameter vector based on the population distribution stored in the
#' \code{ddm} object.
#' @param time_parameters Numeric vector of for setting the time range
#' and time step for simulation (using inverse transform method).
#' @param debug Logical. If \code{TRUE}, print debugging output during
#' simulation. Defaults to \code{FALSE}.
#'
#' @return A data frame with simulated trial data, including columns
#' like \code{subject}, \code{trial}, \code{choice}, and \code{rt}.
#'
#' @details
#' This method simulates data from a design-based Decision Diffusion
#' model (DDM). You can simulate multiple subjects, override default
#' parameters, and use conduct inverse sampling with a fine time step.
#' You may turn on debugging output when something goes wrong.
#'
#' @examples
#' if (requireNamespace("ggdmcModel", quietly = TRUE)) {
#'   BuildModel <- getFromNamespace("BuildModel", "ggdmcModel")
#'   model <- ggdmcModel::BuildModel(
#'     p_map = list(
#'       a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
#'       t0 = "1", st0 = "1", s = "1", precision = "1"
#'     ),
#'     match_map = list(M = list(s1 = "r1", s2 = "r2")),
#'     factors = list(S = c("s1", "s2")),
#'     constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
#'     accumulators = c("r1", "r2"),
#'     type = "fastdm"
#'   )
#' 
#' # Simulate one subject ------------------
#' sub_model <- setDDM(model)
#' p_vector <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = .38)
#' dat <- simulate(sub_model,
#'   nsim = 256, parameter_vector = p_vector,
#'   n_subject = 1
#' )
#' # Simulate multiple subjects ------------------
#' # 1. First build a population distribution to sample many p_vector's
#'
#' pop_mean <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = 0.38)
#' pop_scale <- c(a = 0.05, sz = 0.01, t0 = 0.02, v = .5, z = 0.01)
#' pop_dist <- ggdmcPrior::BuildPrior(
#'   p0    = pop_mean,
#'   p1    = pop_scale,
#'   lower = c(0, 0, 0, -10, 0),
#'   upper = rep(NA, model@npar),
#'   dists = rep("tnorm", model@npar),
#'   log_p = rep(FALSE, model@npar)
#' )
#'
#' # 2. Enter the population distribution to a DDM model
#' pop_model <- setDDM(model, population_distribution = pop_dist)
#'
#' # 3. simulate method will use the distribution to sample new p_vector's
#' # Note: You must enter sensible pop_mean, pop_scale and distribution.
#' hdat <- simulate(pop_model, nsim = 128, n_subject = 32)
#' }
#'
#' @export
#' @method simulate ddm
setMethod(
  "simulate", signature(object = "ddm"),
  function(object, nsim = 4L, seed = NULL, n_subject = 3L,
           parameter_vector = NULL, time_parameters = c(t_min = -1, t_max = 1, dt = 0.001),
           debug = FALSE) {
    #--- Checking if required arguments are provided ---
    if (is.null(object@population_distribution) && is.null(parameter_vector)) {
      stop("Neither population_distribution (a slot in the 'lba' class) nor parameter_vector was found")
    }

    if (!is.null(object@population_distribution) && !is.null(parameter_vector)) {
      stop("You have provided both the population_distribution and the parameter_vector. Which one do you want me to use?")
    }

    ncell <- length(object@model@cell_names)

    if (ncell == 0) {
      stop("Number of cells (ncell) cannot be zero")
    }

    # --- Parameter Preparation ---
    if (!is.null(parameter_vector)) {
      if (is.null(names(parameter_vector))) {
        stop("The parameter_vector must carry name attribute to ")
      }
      name_sorted_p_vector <- parameter_vector[sort(names(parameter_vector))]
      param_matrix <- t(sapply(seq_len(n_subject), function(i) {
        matrix(name_sorted_p_vector, nrow = 1, ncol = object@model@npar)
      }))
    } else {
      # (!is.null(object@population_distribution)) {
      param_matrix <- .prepare_parameter_matrix(
        model = object@model,
        n_subject = n_subject,
        pop_dist = object@population_distribution,
        seed = seed
      )

      # --- Parameter Validation ---
      param_matrix <- .validate_parameters(
        parameter_matrix = param_matrix,
        pop_dist = object@population_distribution,
        rt_model = object
      )
    }

    # --- Simulation ---
    ## Calculate trials per cell with warning if uneven division
    # (nsim %% ncell == 0)
    if (nsim %% ncell == 0) {
      n_trial_per_cell <- nsim / ncell
    } else {
      n_trial_per_cell <- ceiling(nsim / ncell)
      old_nsim <- nsim
      nsim <- n_trial_per_cell * ncell

      message("Warning: n_trial (", old_nsim, ") is not a multiple of the number of cells (", ncell, "). I modified it to ", nsim, ". n_trial per cell = ", n_trial_per_cell, ".\n")
    }

    message("\n[n_trial per condition, n_trial]: [", n_trial_per_cell, ", ", nsim, "]")

    # cat("Going to simulate_lba_trials")
    # print(param_matrix[1, ])

    #--- Seed Handling ---
    # Generate seeds based on number of cores
    if (is.null(seed)) {
      main_seed <- as.integer(Sys.time()) %% 1000000L
    } else {
      main_seed <- seed
    }

    set.seed(main_seed)
    seeds <- sample.int(1e6, n_subject)

    # Print header information
    message("Simulation settings:")
    message("---------------------")
    message("Main seed: ", main_seed)
    message("Number of subjects: ", n_subject)

    # Print seeds in a more readable format
    if (n_subject <= 5) {
      message("\nSeeds for each subject:")
      for (i in seq_len(n_subject)) {
        message("  Subject ", i, ": ", seeds[i])
      }
    } else {
      message("\nSeeds for the first 5 subjects:")
      for (i in seq_len(5)) {
        message("  Subject ", i, ": ", seeds[i])
      }
    }

    t0 <- Sys.time()
    results_list <- lapply(seq_len(n_subject), function(i) {
      trials <- .simulate_ddm_trials(
        rt_model_r = object,
        parameters_r = param_matrix[i, ],
        time_parameters_r = time_parameters,
        n_trial = as.integer(nsim),
        seed = seeds[i],
        debug = debug
      )

      trials$s <- i
      trials
    })
    t1 <- Sys.time()


    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    msg_str <- paste0("Processing time ", "(process model): ")
    message(msg_str, round(proc_time, 3), " secs.")


    # --- Combine Results ---
    results <- do.call("rbind", results_list)
    results$s <- factor(results$s)

    out <- .separate_condition_column(
      results, object@model@factors,
      object@model@accumulators,
      object@model@match_map
    )

    # --- Add Attributes ---
    colnames(param_matrix) <- object@model@pnames
    if (nrow(param_matrix) == 1) {
      attr(out, "parameters") <- param_matrix[1, ]
    } else {
      attr(out, "parameters") <- param_matrix
    }

    attr(out, "main_seed") <- main_seed
    attr(out, "seeds") <- seeds
    return(out)
  }
)
