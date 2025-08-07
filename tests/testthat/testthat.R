Sys.setenv("R_TESTS" = "")
## Workaround for the error,
## "cannot open file 'startup.Rs': No such file or directory" in Windows 10

library(testthat)
library(lbaModel)
library(ggdmc)
library(rtdists)
cat("\nRunning testthat in the directory: ")
cat(getwd(), "\n")


cat("\n========================== Group 0 tests ==========================\n\n")
# test_file(path = "Group0/0_dlba.R")
# test_file(path = "Group0/1_fptpdf_cdf.R")
# test_file(path = "Group0/2_node1.R")
# test_file(path = "Group0/3_node1_3acc.R")
# test_file(path = "Group0/4_rlba_norm.R")
# test_file(path = "Group0/5_node1_3acc_range.R")
# test_file(path = "Group0/6_dlba_all.R")
# test_file(path = "Group0/7_different_time_resolution.R")

cat("\n========================== Group 1 tests ==========================\n\n")
test_file(path = "Group1/0_simulate_error.r")
