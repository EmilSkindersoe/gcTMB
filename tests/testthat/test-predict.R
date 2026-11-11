test_that("predict equals the mean function at the estimates", {
  for (m in models){
    d <- simgauss(m)
    f <- gcTMB(y ~ day, data=d, model=m)
    expect_equal(predict(f), gcmean(m, d$day, bcoef(f)), tolerance=1e-10, label=m)
  }
})

test_that("standard errors equal the delta method", {
  for (m in models){
    d <- simgauss(m)
    f <- gcTMB(y ~ day, data=d, model=m)
    nd <- data.frame(day=c(0, 5, 10, 15, 30))
    p <- predict(f, newdata=nd, se.fit=TRUE)
    b <- unname(bcoef(f))
    se <- sapply(nd$day, function(t){
      g <- numDeriv::grad(function(bb) gcmean(m, t, bb), b)
      sqrt(drop(t(g) %*% vcov(f)[1:3, 1:3] %*% g))
    })
    expect_equal(p$se.fit, se, tolerance=1e-5, label=m)
    # the same as predict.nls with the delta method
    expect_equal(p$fit, unname(predict(nlsfit(m, d), newdata=nd)), tolerance=1e-6)
  }
})

test_that("confidence and prediction intervals of gaussian models", {
  d <- simgauss("logistic")
  f <- gcTMB(y ~ day, data=d)
  nd <- data.frame(day=c(2, 10, 20))
  ci <- predict(f, newdata=nd, interval="confidence")
  pi <- predict(f, newdata=nd, interval="prediction")
  z <- qt(0.975, df.residual(f))
  expect_equal(ci$lwr, ci$fit - z*ci$se.fit)
  expect_equal(pi$upr, pi$fit + z*sqrt(ci$se.fit^2 + residpar(f)^2))
  expect_true(all(pi$lwr < ci$lwr & pi$upr > ci$upr))
})

test_that("prediction intervals for fitted groups and for new groups", {
  d <- simre()
  f <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=d)
  nd <- data.frame(day=c(20, 20), id=c("1", "new"))
  p <- predict(f, newdata=nd, interval="prediction")
  vc <- gcTMB:::recov(f)[[1]]
  # the new group adds the variance of b1 times the derivative of the mean w.r.t. b1, F(t)
  Fz <- predict(f, newdata=nd[2, ], random=FALSE) / bcoef(f)[[1]]
  se <- predict(f, newdata=nd[2, ], se.fit=TRUE)$se.fit
  expect_equal(p$upr[2] - p$fit[2], qnorm(0.975)*sqrt(se^2 + residpar(f)^2 + vc$sd^2*Fz^2), tolerance=1e-6)
  expect_gt(p$upr[2] - p$lwr[2], p$upr[1] - p$lwr[1])
  expect_equal(predict(f, newdata=nd[2, ]), predict(f, newdata=nd[2, ], random=FALSE))
})

test_that("confidence intervals stay in the range of the mean", {
  set.seed(9)
  d <- data.frame(day=rep(days, each=2), n=10)
  d$k <- rbinom(nrow(d), d$n, gcmean("logistic", d$day, c(0.9, 10, 2.5)))
  f <- gcTMB(cbind(k, n - k) ~ day, data=d, family="binomial")
  ci <- predict(f, newdata=data.frame(day=c(0, 10, 24)), interval="confidence")
  expect_true(all(ci$lwr > 0 & ci$upr < 1))
  # parameter draws with b1 > 1 give no binomial observations, with a warning
  pi <- suppressWarnings(predict(f, newdata=data.frame(day=c(0, 10, 24), k=0, n=10), interval="prediction", seed=1))
  expect_true(all(pi$lwr >= 0 & pi$upr <= 1))
})

test_that("predicted parameters are the coefficients, with their standard errors", {
  d <- rbind(cbind(simgauss("logistic", seed=1), trt="a"),
             cbind(simgauss("logistic", seed=2, b=c(60, 12, 2.5)), trt="b"))
  d$trt <- factor(d$trt)
  f <- gcTMB(y ~ day, fform=b1 + b2 ~ trt, data=d)
  nd <- data.frame(trt=c("a", "b"))
  pp <- predict(f, newdata=nd, type="parameters")
  b <- bcoef(f)
  expect_equal(colnames(pp), c("b1", "b2", "b3"))
  expect_equal(unname(pp[2, ]), unname(c(b[1] + b[2], b[3] + b[4], b[5])))
  ps <- predict(f, newdata=nd, type="parameters", se.fit=TRUE)
  V <- vcov(f)
  expect_equal(ps$se.fit[ps$parameter == "b1"], sqrt(c(V[1, 1], V[1, 1] + V[2, 2] + 2*V[1, 2])), tolerance=1e-8)
  # the logistic b2 is the time of timeto at 50%
  expect_equal(unname(pp[, "b2"]), unname(coef(timeto(f, newdata=nd))))
})

test_that("bootstrap intervals of predictions", {
  d <- simgauss("logistic")
  f <- gcTMB(y ~ day, data=d)
  bt <- bootstrap(f, nsim=100, seed=1, progress=FALSE)
  nd <- data.frame(day=c(5, 10, 20))
  pb <- predict(f, newdata=nd, interval="confidence", boot=bt)
  pd <- predict(f, newdata=nd, interval="confidence")
  expect_equal(pb$fit, pd$fit)
  expect_equal(pb$se.fit, pd$se.fit, tolerance=0.3)
})
