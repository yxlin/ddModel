q(save = "no")
cat("\n\n-------------------- rfastdm --------------------")
rm(list = ls())
pkg <- c("ddModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

# /* p_vector, passing also RT and precision
#      * 0  1   2   3   4   5   6  7   8   9
#      * a  v  zr   d  sz  sv  t0 st0  s   precision
#      */
# RT <- seq(0.001, 1.2, 0.001)
# RT <- 0.550
## no var ------------------
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, szr = 0, sv = 0.0, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
set.seed(123)
time_parameters <- c(-1, 1, 0.01)
res0 <- rfastdm(n = 1, parameters_r = p_vector, time_parameters_r = time_parameters, debug = T)
res0


# set.seed(123)
# res0 <- rfastdm(p_vector, n = 32)
# res0
