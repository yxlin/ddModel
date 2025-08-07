#' The Decision Diffusion Model
#'
#' \pkg{ddModel} provides functions for computing the density, distribution,
#' and random generation of the Decision Diffusion model (DDM), a widely used
#' cognitive model for analysing choice and response time data. The package
#' allows model specification, including the ability to fix, constrain, or vary
#' parameters across experimental conditions. While it does not include a
#' built-in optimiser, it supports likelihood evaluation and can be integrated
#' with external tools for parameter estimation. Functions for simulating
#' synthetic datasets aer also provided. This package is intended for
#' researchers modelling speeded decision-making in behavioural and cognitive
#' experiments. For more information, see
#' Voss, Rothermund, and Voss (2004) <doi:10.3758/BF03196893>,
#' Voss and Voss (2007) <doi:10.3758/BF03192967>, and
#' Ratcliff and McKoon (2008) <doi:10.1162/neco.2008.12-06-420>.
#'
#' @keywords package
#'
#' @name ddModel
#' @keywords internal
#' @author  Yi-Shin Lin <yishinlin001@gmail.com> \cr
#' @importFrom Rcpp evalCpp
#' @useDynLib ddModel
"_PACKAGE"
NULL

#' Density, Distribution and Random Generation of the Decision Diffusion Model
#'
#' A set of functions implementing the Decision Diffusion Model (DDM),
#' providing: probability density (\code{dfastdm}), cumulative distribution
#' (\code{pfastdm}), and random choices and response times
#' generation (\code{rfastdm}).
#'
#' @name fastdm
#' @aliases dfastdm pfastdm rfastdm
#'
#' @param rt_r A numeric vector of response times (in seconds)
#' @param parameters_r A numeric vector of model parameters:
#' \describe{
#'   \item{a}{Boundary separation.}
#'   \item{z}{Starting point (bias), must lie between 0 and\code{a}.}
#'   \item{v}{Drift rate.}
#'   \item{t0}{Non-decision time (must be \eqn{\ge 0}).}
#'   \item{d}{Difference in motor execution time between boundaries (e.g.,
#'            left vs. right response). When \code{d > 0}, responses at the
#'            upper boundary are delayed by \code{d} seconds. Symmetric 
#'            execution is assumed when \code{d = 0} (Voss and Voss, 2007.)}
#'   \item{sz}{Inter-trial variability in starting point (optional)}
#'   \item{sv}{Inter-trial variability in drift rate (optional)}
#'   \item{st0}{Inter-trial variability in non-decision time (optional)}
#' }
#' 
#' @param is_lower Logical; if \code{TRUE}, compute probability for the lower
#' boundary. Defaults to \code{TRUE}. Note: under standard DDM assumptions,
#' a response at the upper boundary indicates a 'correct' choice.
#' @param n Integer; number of random samples to generate (used in
#' \code{rfastdm}).
#' @param time_parameters_r Optional numeric vector specifying time resolution
#' for simulation (used in \code{rfastdm}, which employs inverse transform
#' sampling).
#' @param debug Logical; if \code{TRUE}, return debugging information.
#'
#' @return
#' \describe{
#'   \item{dfastdm}{Numeric vector of probability densities}
#'   \item{pfastdm}{Numeric vector of cumulative probabilities}
#'   \item{rfastdm}{DataFrame with columns: \code{rt} (response time),
#'                 \code{response} (0 = lower/incorrect, 1 = upper/correct).}
#' }
#'
#' @details
#' These functions implement the DDM, which describes a one-dimensional
#' evidence accumulation process between two absorbing boundaries, representing
#' competing response alternatives. The implementation includes an extended
#' version following the \code{fast-dm} framework, allowing for inter-trial
#' variability in key parameters (e.g., drift rate, starting point,
#' non-decision time).
#'
#' The density and distribution functions are based on the \code{g_minus} and
#' \code{g_plus} formulations (in fastdm source code), which were described
#' in Equation A1–A4 in Voss, Rothermund, and Voss  (2004). These equations
#' described analytic solutions to Ratcliff’s (1978) diffusion model. The
#' original C implementation was part of the \code{fast-dm 30.2}.
#' This version was rewritten in modern C++ and adapted for use in Differential
#' Evolution Markov Chain Monte Carlo (DE-MCMC) estimation.
#'
#' @references
#' Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the parameters of
#' the diffusion model: An empirical validation. \emph{Memory & Cognition},
#' \bold{32}(7), 1206–1220.
#'
#' Voss, A., & Voss, J. (2007). Fast-dm: A free program for efficient diffusion
#' model analysis. \emph{Behavior Research Methods}, \bold{39}, 767–775.
#' \doi{10.3758/BF03192967}
#'
#' @examples
#' # Density function
#' RT <- 0.550
#' unsorted_p_vector <- c(
#'     a = 1, v = 1.5, z = 0.5 * 1, d = 0, sz = 0, sv = 0, t0 = 0.15, st0 = 0,
#'     s = 1, precision = 3
#' )
#' p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
#'
#' result0 <- dfastdm(RT, p_vector)
#'
#'
#' # Cumulative distribution function
#' RT <- seq(0.1, 1.2, 0.01)
#' unsorted_p_vector <- c(
#'     a = 1, v = 1.5, zr = 0.5 * 1, d = 0, sz = 0.0,
#'     sv = 0, t0 = 0.15, st0 = 0, s = 1, precision = 3
#' )
#' p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
#'
#' result1 <- pfastdm(RT, p_vector)
#'
#' sz <- 0.05
#' sv <- 0.01
#' st0 <- 0.001
#'
#' unsorted_p_vector <- c(
#'     a = 1, v = 1.5, zr = 0.5 * 1, d = 0, sz = sz, sv = sv, t0 = 0.15,
#'     st0 = st0,
#'     s = 1, precision = 3
#' )
#' p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
#' result2 <- pfastdm(RT, p_vector, is_lower = TRUE)
#'
#' # Random generation
#' unsorted_p_vector <- c(
#'     a = 1, v = 1.5, zr = 0.5 * 1, d = 0, szr = 0, sv = 0.0, t0 = 0.15,
#'     st0 = 0, s = 1, precision = 3
#' )
#' p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
#' set.seed(123)
#' time_parameters <- c(-1, 1, 0.01)
#' result3 <- rfastdm(
#'     n = 1, parameters_r = p_vector,
#'     time_parameters_r = time_parameters, debug = TRUE
#' )
#'
#' # Debugging information
#' # st0 = 0 < 1.50017e-08. sv = 0 < 1e-05. sz = 0 < 1e-05. Selecting f_plain.
#' #
#' # Initial search range: t_min = -1, t_max = 1
#' # Point 0: t = -1, F = 0.00139754
#' # Point 1: t = -0.99, F = 0.00148506
#' # Point 2: t = -0.98, F = 0.0015778
#' # Point 3: t = -0.97, F = 0.00167698
#' # Point 4: t = -0.96, F = 0.00178253
#' # Point 196: t = 0.96, F = 0.99201
#' # Point 197: t = 0.97, F = 0.992483
#' # Point 198: t = 0.98, F = 0.992927
#' # Point 199: t = 0.99, F = 0.993343
#' # Point 200: t = 1, F = 0.993735
NULL
