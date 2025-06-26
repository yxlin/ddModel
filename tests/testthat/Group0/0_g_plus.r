# q(save = "no")
cat("\n\n-------------------- g_plus --------------------")
rm(list = ls())
pkg <- c("ddModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

# /* p_vector, passing also RT and precision
#      * 0  1   2   3   4   5   6  7   8   9
#      * a  v  zr   d  sz  sv  t0 st0  s   precision
#      */

## no var ------------------

RT <- 0.550

# x <- seq(0, 3, length.out = 1e2)
unsorted_p_vector <- c(
    a = 1, v = 1.5, z = 0.5 * 1, d = 0, sz = 0, sv = 0, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]

res0 <- dfastdm(RT, p_vector)
res1 <- rtdiststest::ddiffusion(RT, a = 1, v = 1.5, t0 = .15, response = "lower")
cat("Lower boundary - Incorrect: ", all.equal(res0, res1))


RT <- seq(0.001, 1.2, 0.001)

res0 <- dfastdm(RT, p_vector, is_lower = FALSE)
res1 <- rtdiststest::ddiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = 0, response = "upper"
)
cat("Upper boundary - Correct: ", all.equal(res0, res1))

## sv, sz = 0 ------------------
RT <- seq(0.001, 1.2, 0.001)

sz <- 0
sv <- 0.1
unsorted_p_vector <- c(
    a = 1, v = 1.5, z = 0.5 * 1, d = 0, sz = sz, sv = sv, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]

res0 <- dfastdm(RT, p_vector, is_lower = FALSE)
res1 <- rtdiststest::ddiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = sz, sv = sv,
    st0 = 0, response = "upper"
)

cat("sz = 0, sv != 0, Upper boundary - Correct: ", all.equal(res0, res1))


## szr ------------------
x <- seq(0, 3, length.out = 1e2)
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, szr = 0.5, sv = 0, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]

res0 <- dfastdm(RT, p_vector)

res1 <- rtdiststest::ddiffusion(RT, a = 1, v = 1.5, t0 = .15, sz = 0.5, response = "lower")

all.equal(res0, res1)

## sv ------------------
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, szr = 0.5, sv = 0.1, t0 = 0.15, st0 = 0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
res0 <- dfastdm(RT, p_vector)

res1 <- rtdiststest::ddiffusion(RT, a = 1, v = 1.5, t0 = .15, sz = 0.5, sv = 0.1, response = "lower")


all.equal(res0, res1)

## st0 ------------------
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, szr = 0.5, sv = 0.1, t0 = 0.15, st0 = 0.1,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]
res0 <- dfastdm(RT, p_vector)
res1 <- rtdiststest::ddiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = 0.5, sv = 0.1,
    st0 = 0.1, response = "lower"
)
all.equal(res0, res1)


sz <- 0.1
sv <- 0.1
st0 <- 0.1
unsorted_p_vector <- c(
    a = 1, v = 1.5, zr = 0.5 * 1, d = 0, szr = sz, sv = sv, t0 = 0.15, st0 = st0,
    s = 1, precision = 3
)
p_vector <- unsorted_p_vector[sort(names(unsorted_p_vector))]

res0 <- dfastdm(RT, p_vector)
res1 <- rtdiststest::ddiffusion(RT,
    a = 1, v = 1.5, t0 = .15, sz = sz, sv = sv,
    st0 = st0, response = "lower"
)

cat("Lower boundary - Incorrect: ", all.equal(res0, res1), "\n")
