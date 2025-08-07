q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ddModel", "ggdmcPrior", "ggdmcPhi")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ddModel/tests/testthat/Group1/"

p_map <- list(
    a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1", t0 = "1",
    st0 = "1", s = "1"
)
factors <- list(S = c("s1", "s2"), F = c("f1", "f2"))
accumulators <- c("r1", "r2")
match_map <- list(M = list(s1 = "r1", s2 = "r2"))

parameter_M <- ggdmcModel::bind_condition2parameters_r(p_map, factors)
parameter_M

cell_and_factor_names <- ggdmcModel::build_cell_names_r(p_map, factors, accumulators)
cell_and_factor_names[[1]]
cell_and_factor_names[[2]]

model_boolean <- ggdmcModel::build_model_boolean_r(p_map, factors, accumulators, match_map)


model_r1 <- model_boolean[, , 1]
colnames(model_r1) <- parameter_M
rownames(model_r1) <- cell_and_factor_names[[1]]
model_r1
