Sys.setenv("R_TESTS" = "")
## Workaround for the error,
## "cannot open file 'startup.Rs': No such file or directory" in Windows 10

library(testthat)
library(ddModel)
cat("\nRunning testthat in the directory: ")
cat(getwd(), "\n")


