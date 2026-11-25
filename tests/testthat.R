library(testthat)
library(gcTMB)

# GCTMB_TEST_REPORTER=location (set for one row of the R-CMD-check workflow) reports each test as it
# starts, on stderr (unbuffered), so that the output shows the test during which R crashes
if (Sys.getenv("GCTMB_TEST_REPORTER") == "location"){
  test_check("gcTMB", reporter=LocationReporter$new(file=stderr()))
} else {
  test_check("gcTMB")
}
