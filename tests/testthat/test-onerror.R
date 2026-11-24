test_that("control onerror='warning' returns a fit whose optimiser fails", {
  d <- FabaHeights25
  expect_error(gcTMBcontrol(onerror="ignore"))
  expect_equal(gcTMBcontrol()$onerror, "stop")

  local_mocked_bindings(nlminb=function(...) stop("optimiser broke"))
  expect_error(gcTMB(height ~ DAS, data=d), "optimiser broke")
  expect_warning(f <- gcTMB(height ~ DAS, data=d, control=list(onerror="warning")),
                 "the fit failed: optimiser broke; the object is returned at the starting values")
  expect_s3_class(f, "gcTMB")
  expect_match(f$fiterror, "optimiser broke")
  # at the starting values, without standard errors
  expect_equal(unname(f$estimates[seq_along(f$start[!f$fix])]), unname(f$start[!f$fix]))
  expect_true(all(is.na(f$ssdrl[, 2])))
  cc <- converged(f)
  expect_false(cc)
  expect_match(attr(cc, "reasons")[1], "the fit failed: optimiser broke")
  expect_output(print(f), "The fit failed: optimiser broke")
  expect_output(print(summary(f)), "The fit failed: optimiser broke")

  # with random effects
  expect_warning(fr <- gcTMB(height ~ DAS, rform=b1 ~ 1 | plant_ID, data=d, control=list(onerror="warning")),
                 "the fit failed")
  expect_false(converged(fr))
})

test_that("control onerror='warning' returns a fit whose sdreport fails, without standard errors", {
  d <- FabaHeights25
  f0 <- gcTMB(height ~ DAS, data=d)
  local_mocked_bindings(sdreport=function(...) stop("no covariance"))
  expect_error(gcTMB(height ~ DAS, data=d), "no covariance")
  expect_warning(f <- gcTMB(height ~ DAS, data=d, control=list(onerror="warning")),
                 "returned without standard errors")
  # the estimates of the successful optimisation
  expect_equal(coef(f), coef(f0), tolerance=1e-8)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f0)))
  expect_true(all(is.na(f$ssdrl[, 2])))
  expect_false(converged(f))
})

test_that("summary repeats the warning of a fit that did not converge", {
  expect_warning(f <- gcTMB(height ~ DAS, rform=b1 ~ 1 | plant_ID, data=FabaHeights25, control=list(iter.max=2)),
                 "the optimiser did not converge")
  expect_output(print(summary(f)), paste("Warning: the optimiser did not converge:", f$opt$message), fixed=TRUE)
  # a converged fit has no such warning
  s <- capture.output(print(summary(gcTMB(height ~ DAS, data=FabaHeights25))))
  expect_false(any(grepl("did not converge", s)))
})
