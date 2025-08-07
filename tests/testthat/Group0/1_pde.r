q(save = "no")
cat("\n\n-------------------- pde --------------------")
rm(list = ls())
pkg <- c("ddModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

# /* p_vector, passing also RT and precision

#      * 0  1   2   3   4   5   6  7   8   9
#      * a  v  zr   d  sz  sv  t0 st0  s   precision
#      */


# RT <- c(0.1, 0.550)
# RT <- seq(0.1, 1.2, 0.1)
RT <- seq(0.1, 1.2, 0.01)
# RT <- 0.11
# RT <- RT[2:10]
# RT <- RT[6:7]

# no var ----------------
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, sz = 0.0, sv = 0, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]

res0 <- pfastdm(RT, p_vector)

# cat("\n")
res1 <- rtdiststest::pdiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = 0.0,
    sv = 0.0, st0 = 0, s = 1, precision = 3, response = "lower"
)
all.equal(res0, res1)


# sz ----------------
unsorted_p_vector <- c(
    a = 1, v = 1.5, z = 0.5 * 1, d = 0, sz = 0.1, sv = 0, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
p_vector
res0 <- pfastdm(RT, p_vector)
res1 <- rtdiststest::pdiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = 0.1,
    sv = 0.0, st0 = 0, s = 1, precision = 3, response = "lower"
)

all.equal(res0, res1)
# plot(RT, res0)

# sv ----------------
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, sz = 0.1, sv = 0.1, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
res0 <- pfastdm(RT, p_vector)
res1 <- rtdiststest::pdiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = 0.1,
    sv = 0.1, st0 = 0, s = 1, precision = 3, response = "lower"
)
all.equal(res0, res1)
# plot(RT, res0)

# st0 ----------------
sz <- 0.05
sv <- 0.01
st0 <- 0.001

unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, sz = sz, sv = sv, t0 = 0.15,
    st0 = st0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
res0 <- pfastdm(RT, p_vector, is_lower = TRUE, debug = TRUE)
cat("\n\n rtdists\n")

res1 <- rtdiststest::pdiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = sz,
    sv = sv, st0 = st0, s = 1, precision = 3, response = "lower"
)

all.equal(res0, res1)
