# q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ddModel", "ggdmcPrior", "ggdmcPhi")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ddModel/tests/testthat/Group1/"

# p_map <- list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1", s = "1")
# factors <- list(S = c("s1", "s2"))
# parameter_M <- ggdmcModel::bind_condition2parameters_r(p_map, factors)
# parameter_M

# res <- ggdmcModel::is_core_parameter_x_condition(p_map, factors)
# names(res) <- sort(names(p_map))
# res


p_map <- list(a = "1", v = "S", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1")
factors <- list(S = c("s1", "s2"))
parameter_M <- ggdmcModel::bind_condition2parameters_r(p_map, factors)
parameter_M

res0 <- ggdmcModel::is_core_parameter_x_condition(p_map, factors)
res1 <- ggdmcModel::is_parameter_x_condition(p_map, factors)
names(res0) <- sort(names(p_map))
names(res1) <- parameter_M
res0
res1

# LBA st0 fix -------------------
# p_map <- list(A = "1", B = "1", t0 = "D", mean_v = "M", sd_v = "1", st0 = "1")
# factors <- list(S = c("red", "blue"), D = c("easy", "hard"))
# parameter_M <- ggdmcModel::bind_condition2parameters_r(p_map, factors)
# res0 <- ggdmcModel::is_core_parameter_x_condition(p_map, factors)
# res1 <- ggdmcModel::is_parameter_x_condition(p_map, factors)
# names(res0) <- sort(names(p_map))
# names(res1) <- parameter_M
# res0
# res1

# res <- ggdmcModel::split_parameter_x_condition(parameter_M)
# unlist(res)

# res <- ggdmcModel::get_stimulus_level_r(p_map, factors, accumulators)
# res

# res <- ggdmcModel::get_factor_cells_r(p_map, factors, accumulators)
# unlist(res)


# cell_and_factor_names <- ggdmcModel::build_cell_names_r(p_map, factors, accumulators)
# cell_and_factor_names[[1]]
# cell_and_factor_names[[2]]

# parameter_x_condition_names <- ggdmcModel::bind_condition2parameters_r(p_map, factors)
# parameter_x_condition_names

# p_map <- list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1", s = "1")
# factors <- list(S = c("s1", "s2"), F = c("f1", "f2"))
# accumulators <- c("r1", "r2")
# match_map <- list(M = list(s1 = "r1", s2 = "r2"))

# model_boolean <- ggdmcModel::build_model_boolean_r(p_map, factors, accumulators, match_map)

# hyper_model@model_boolean


# hyper_model <- ggdmcModel::BuildModel(
#     p_map = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1"),
#     match_map = list(M = list(s1 = "r1", s2 = "r2")),
#     factors = list(S = c("s1", "s2"), F = c("f1", "f2")),
#     constants = c(d = 1, st0 = 0),
#     accumulators = c("r1", "r2"),
#     type = "fastdm"
# )
