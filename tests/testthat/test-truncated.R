test_that("truncated curves are the truncated distribution functions, zero up to the threshold", {
  t <- seq(0, 30, 0.5)
  for (mod in 1:2){
    F <- gcTMB:::gcmodel(mod)$F
    f <- gcTMB:::gcmeanR(mod, t, 10, 12, 3, 6, 1, trunc=TRUE)
    naive <- ifelse(t > 6, 10 * (F((t - 12)/3) - F((6 - 12)/3)) / (1 - F((6 - 12)/3)), 0)
    expect_equal(f, naive, tolerance=1e-12)
    expect_true(all(f[t <= 6] == 0))
    # without truncation in the limit of a threshold far below b2
    expect_equal(gcTMB:::gcmeanR(mod, t, 10, 12, 3, -200, 1, trunc=TRUE),
                 gcTMB:::gcmeanR(mod, t, 10, 12, 3, 0, 1), tolerance=1e-10)
    # thresholds far beyond b2, where 1 - F(z4) underflows
    expect_equal(gcTMB:::gcmeanR(mod, 81, 10, 12, 1, 80, 1, trunc=TRUE), 10 * (1 - exp(-1)), tolerance=1e-6)
    # times to reach a level are the inverse of the curve
    q <- c(0.1, 0.5, 0.9)
    tt <- gcTMB:::gctimeR(mod, q, 12, 3, 6, trunc=TRUE)
    expect_equal(gcTMB:::gcmeanR(mod, tt, 1, 12, 3, 6, 1, trunc=TRUE), q, tolerance=1e-12)
  }
})

simtrunc <- function(mod, seed=2){
  set.seed(seed)
  d <- data.frame(day=rep(0:30, each=4))
  d$y <- gcTMB:::gcmeanR(mod, d$day, 50, 12, 3, 8, 1, trunc=TRUE) + rnorm(nrow(d), 0, 1)
  d
}

test_that("truncated logistic and gompertz curves recover the threshold", {
  for (model in c("logistic", "gompertz")){
    d <- simtrunc(match(model, c("logistic", "gompertz")))
    fit <- gcTMB(y ~ day, data=d, model=model, fix=c(b4=FALSE))
    expect_equal(fit$opt$convergence, 0)
    cf <- coef(fit)
    expect_equal(unname(cf[["b4:(Intercept)"]]), 8, tolerance=0.1)
    expect_equal(unname(cf[["b1:(Intercept)"]]), 50, tolerance=0.05)
    # better than the curve without truncation
    expect_lt(AIC(fit), AIC(gcTMB(y ~ day, truncated=FALSE, data=d, model=model)))
  }
})

test_that("truncated fits: gradient, predictions and times", {
  d <- simtrunc(1)
  fit <- gcTMB(y ~ day, data=d, fix=c(b4=FALSE))
  expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-6)
  expect_equal(as.vector(fit$obj$gr(fit$opt$par)), numDeriv::grad(fit$obj$fn, fit$opt$par), tolerance=1e-5)
  nd <- data.frame(day=c(5, 9, 15))
  p <- predict(fit, newdata=nd, se.fit=TRUE)
  expect_equal(unname(unlist(p$fit)), unname(predict(fit, newdata=nd)), tolerance=1e-10)
  expect_equal(unlist(p$fit)[[1]], 0)
  # the curve reaches the levels at the times of timeto
  tt <- timeto(fit, level=c(10, 50, 90))
  b1 <- coef(fit)[["b1:(Intercept)"]]
  expect_equal(unname(predict(fit, newdata=data.frame(day=coef(tt)))), c(0.1, 0.5, 0.9) * b1, tolerance=1e-8)
  expect_true(all(is.finite(sqrt(diag(vcov(tt))))))
  expect_output(print(fit), "truncated at the threshold b4")
})

test_that("zero counts at a mean of zero before the threshold", {
  # the maximum is at a kink, at the observation time 65 with only zero counts
  expect_no_warning(fit <- gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE)))
  expect_equal(fit$opt$convergence, 0)
  expect_match(fit$opt$message, "kink")
  expect_equal(coef(fit)[["b4:(Intercept)"]], 65)
  # the maximum with b4 held at 65, the same with log links
  fit65 <- gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson", start=c(b4=65), fix=c(b4=TRUE))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(fit65)), tolerance=1e-8)
  fitl <- gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE), link=c(b1="log", b4="log"))
  expect_equal(as.numeric(logLik(fitl)), as.numeric(logLik(fit)), tolerance=1e-8)
  expect_true(is.finite(logLik(fit)))
  b4 <- coef(fit)[["b4:(Intercept)"]]
  expect_true(all(FabaCounts23$pods[FabaCounts23$DAS <= b4] == 0))
  r <- residuals(fit, type="pearson")
  expect_true(all(is.finite(r)))
  expect_true(all(r[FabaCounts23$DAS <= b4] == 0))
  expect_true(is.finite(summary(fit)$dispersion[["ratio"]]))
  s <- simulate(fit, nsim=2, seed=1)
  expect_true(all(unlist(s)[rep(FabaCounts23$DAS <= b4, 2)] == 0))
})

test_that("b4 is only used by truncated curves", {
  expect_error(gcTMB(height ~ DAS, truncated=FALSE, rform=b4 ~ 1 | plant_ID, data=FabaHeights25), "does not use b4")
  expect_error(gcTMB(height ~ DAS, data=FabaHeights25, truncated=NA), "truncated must be TRUE or FALSE")
  expect_false("b4:(Intercept)" %in% names(coef(gcTMB(height ~ DAS, truncated=FALSE, data=FabaHeights25))))
})

test_that("by default the curves are truncated at b4 = 0, fixed", {
  f <- gcTMB(height ~ DAS, data=FabaHeights25)
  expect_true(f$truncated)
  expect_equal(f$fix, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(f$start[4:5], c(0, 1))
  out <- capture.output(print(f))
  expect_true(any(grepl("truncated at b4 = 0", out)))
  expect_false(any(grepl("^Fixed:", out)))
  expect_equal(colnames(predict(f, type="parameters")), c("b1", "b2", "b3"))
  # the logistic curve starts at zero at t = 0
  b <- bcoef(f)
  expect_equal(gcTMB:::gcmeanR(1, c(0, 20), b[[1]], b[[2]], b[[3]], 0, 1, trunc=TRUE),
               b[[1]]*(plogis((c(0, 20) - b[[2]])/b[[3]]) - plogis(-b[[2]]/b[[3]]))/(1 - plogis(-b[[2]]/b[[3]])))
  # a threshold at zero does not truncate the hill and weibull curves
  for (m in c("hill", "weibull")){
    expect_equal(gcTMB(height ~ DAS, data=FabaHeights25, model=m)$opt$objective,
                 gcTMB(height ~ DAS, data=FabaHeights25, model=m, truncated=FALSE)$opt$objective, tolerance=1e-8)
  }
  # an estimated threshold, and the log link of b4 needs a value
  fb <- gcTMB(height ~ DAS, data=FabaHeights25, fix=c(b4=FALSE))
  expect_equal(fb$fix, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_true(any(grepl("truncated at the threshold b4", capture.output(print(fb)))))
  expect_error(gcTMB(height ~ DAS, data=FabaHeights25, link=c(b4="log")), "outside the range of its log link")
})
