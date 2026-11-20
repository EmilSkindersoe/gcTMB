# The normal and gamma distribution functions and the linear model:
#   normal: b1 Phi((t - b2)/b3)
#   gamma:  b1 P(b3, t (b3 - 1)/b2), shape b3 and scale b2/(b3 - 1), zero up to t = 0
#   linear: b1 t + b6

test_that("the curves are the normal and gamma distribution functions and the line", {
  t <- c(-2, 0, 5, 20, 30, 45, 80)
  tp <- pmax(t, 0)
  expect_equal(gcTMB:::gcmeanR(5, t, 50, 30, 6, 0), 50*pnorm((t - 30)/6))
  expect_equal(gcTMB:::gcmeanR(5, t, 50, 30, 6, 0, 2), 50*pnorm((t - 30)/6)^2)
  expect_equal(gcTMB:::gcmeanR(6, t, 50, 30, 4, 0), ifelse(t > 0, 50*pgamma(tp, 4, scale=10), 0))
  expect_equal(gcTMB:::gcmeanR(6, t, 50, 30, 4, 0, 0.5), ifelse(t > 0, 50*pgamma(tp, 4, scale=10)^0.5, 0))
  expect_equal(gcTMB:::gcmeanR(7, t, 2, 0, 0, 0, f6=5), 2*t + 5)
  # truncated curves
  expect_equal(gcTMB:::gcmeanR(5, t, 50, 30, 6, 10, trunc=TRUE),
               ifelse(t > 10, 50*(pnorm((t - 30)/6) - pnorm(-20/6))/(1 - pnorm(-20/6)), 0))
  # times to reach a level are the inverse of the curves, also truncated and asymmetric
  q <- c(0.1, 0.5, 0.9)
  for (mod in 5:6){
    for (a in c(0.5, 1, 2)){
      for (trunc in c(FALSE, TRUE)){
        tt <- gcTMB:::gctimeR(mod, q, 30, 4, 10, a, trunc=trunc)
        expect_equal(gcTMB:::gcmeanR(mod, tt, 1, 30, 4, 10, a, trunc=trunc), q, tolerance=1e-10)
      }
    }
  }
  expect_equal(gcTMB:::gctimeR(7, (50 - 5)/2, 0, NA, b1=2), (50 - 5)/2)
})

test_that("b2 is the inflection point of the symmetric normal and gamma curves", {
  expect_equal(gcTMB:::maxrateR(5, c(50, 30, 6, 0, 1)), c(tmax=30, rate=50*dnorm(0)/6, fraction=0.5))
  expect_equal(gcTMB:::maxrateR(6, c(50, 30, 4, 0, 1)), c(tmax=30, rate=50*dgamma(30, 4, scale=10), fraction=pgamma(3, 4)))
  # asymmetric curves: the maximum of the slope on a fine grid
  for (mod in 5:6){
    for (a in c(0.5, 2.5)){
      t <- seq(0.01, 150, by=0.001)
      f <- gcTMB:::gcmeanR(mod, t, 50, 30, 4, 0, a)
      r <- diff(f)/0.001
      i <- which.max(r)
      expect_equal(unname(gcTMB:::maxrateR(mod, c(50, 30, 4, 0, a))), c(t[i] + 0.0005, r[i], mean(f[i + 0:1])/50),
                   tolerance=1e-3, label=paste(mod, a))
    }
  }
  expect_equal(gcTMB:::maxrateR(7, c(2, 0, 0, 0, 1, 5)), c(tmax=NA, rate=2, fraction=NA))
})

test_that("normal and gamma fits agree with nls", {
  set.seed(5)
  d <- data.frame(day=rep(seq(2, 60, by=4), 3))
  d$y <- 50*pnorm((d$day - 25)/8) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, data=d, model="normal", truncated=FALSE)
  n <- nls(y ~ b1*pnorm((day - b2)/b3), data=d, start=list(b1=50, b2=25, b3=8))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
  d$y <- 50*pgamma(d$day, 5, scale=25/4) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, data=d, model="gamma")
  n <- nls(y ~ b1*pgamma(day*(b3 - 1)/b2, b3), data=d, start=list(b1=50, b2=25, b3=5))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
  # timeto of the gamma curve, with a profile interval
  b <- bcoef(f)
  tt <- timeto(f, level=50)
  expect_equal(unname(coef(tt)), qgamma(0.5, b[[3]], scale=b[[2]]/(b[[3]] - 1)), tolerance=1e-8)
  ci <- suppressMessages(confint(tt))
  expect_true(ci[1] < coef(tt) && ci[2] > coef(tt))
})

test_that("the linear model is least squares, without truncation and without an asymptote", {
  f <- gcTMB(height ~ DAS, data=FabaHeights25, model="linear")
  l <- lm(height ~ DAS, data=FabaHeights25)
  expect_equal(unname(bcoef(f)), unname(rev(coef(l))), tolerance=1e-6)
  expect_false(f$truncated)
  expect_equal(f$fix, c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_false(any(grepl("truncated", capture.output(print(f)))))
  expect_equal(unname(sqrt(diag(vcov(f))[1:2])), unname(rev(sqrt(diag(vcov(l))))), tolerance=1e-5)
  tt <- timeto(f, level=50, type="absolute")
  expect_equal(unname(coef(tt)), (50 - coef(l)[[1]])/coef(l)[[2]], tolerance=1e-6)
  expect_error(timeto(f, level=50), "no asymptote")
  expect_message(confint(tt), "not available for the linear model")
  mr <- maxrate(f)
  expect_equal(coef(mr)[["rate"]], bcoef(f)[[1]])
  expect_equal(confint(mr, parm="rate", method="wald"), confint(f, parm=1, method="wald"), ignore_attr=TRUE)
})
