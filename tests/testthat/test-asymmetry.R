# The asymmetry b5: Richards (logistic), b1 / (1 + b5 exp(-(t - b2)/b3))^(1/b5); hill as LL.5 of drc,
# b1 / (1 + (b2/t)^b3)^b5; weibull, b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3))^b5; gompertz, a power of the
# survival function, b1 (1 - (1 - exp(-exp(-(t - b2)/b3)))^b5). b5 = 1 gives the symmetric curves.

test_that("the asymmetric curves are the Richards, LL.5, exponentiated weibull and gompertz functions", {
  t <- c(-1, 0, 2, 10, 25, 40, 80)
  tp <- pmax(t, 0)
  a <- 2.5
  expect_equal(gcTMB:::gcmeanR(1, t, 50, 20, 5, 0, a), 50/(1 + a*exp(-(t - 20)/5))^(1/a))
  expect_equal(gcTMB:::gcmeanR(2, t, 50, 20, 5, 0, a), 50*(1 - (1 - exp(-exp(-(t - 20)/5)))^a))
  expect_equal(gcTMB:::gcmeanR(3, t, 50, 20, 4, 0, a), ifelse(t > 0, 50/(1 + (20/tp)^4)^a, 0))
  expect_equal(gcTMB:::gcmeanR(4, t, 50, 20, 3, 0, a), ifelse(t > 0, 50*(1 - exp(-(2/3)*(tp/20)^3))^a, 0))
  # b5 = 1 is the symmetric curve, and the Richards curve tends to the gompertz curve as b5 -> 0
  for (mod in 1:4) expect_equal(gcTMB:::gcmeanR(mod, t, 50, 20, 4, 0, 1), gcTMB:::gcmeanR(mod, t, 50, 20, 4, 0))
  expect_equal(gcTMB:::gcmeanR(1, t, 50, 20, 5, 0, 1e-7), gcTMB:::gcmeanR(2, t, 50, 20, 5, 0, 1), tolerance=1e-5)
})

test_that("times to reach a level are the inverse of the asymmetric curves, also truncated", {
  q <- c(0.05, 0.37, 0.5, 0.95)
  for (mod in 1:4){
    for (a in c(0.3, 1, 3)){
      for (trunc in c(FALSE, TRUE)){
        tt <- gcTMB:::gctimeR(mod, q, 20, 4, 6, a, trunc=trunc)
        expect_equal(gcTMB:::gcmeanR(mod, tt, 1, 20, 4, 6, a, trunc=trunc), q, tolerance=1e-10)
      }
    }
  }
})

test_that("an estimated asymmetry agrees with nls, and timeto with the curve", {
  set.seed(3)
  d <- data.frame(day=rep(seq(2, 60, by=4), 3))
  d$y <- gcTMB:::gcmeanR(1, d$day, 50, 20, 5, 0, 3) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, truncated=FALSE, data=d, fix=c(b5=FALSE))
  expect_equal(f$fix, c(FALSE, FALSE, FALSE, TRUE, FALSE))
  n <- nls(y ~ b1/(1 + b5*exp(-(day - b2)/b3))^(1/b5), data=d, start=list(b1=50, b2=20, b3=5, b5=3))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-4)
  b <- bcoef(f)
  tt <- timeto(f, respLev=c(10, 50, 90))
  expect_equal(unname(coef(tt)), gcTMB:::gctimeR(1, c(0.1, 0.5, 0.9), b[[2]], b[[3]], b5=b[[4]]), tolerance=1e-8)
  # the predictions of the template are the curve
  nd <- data.frame(day=c(5, 20, 50))
  expect_equal(as.vector(predict(f, newdata=nd)), gcTMB:::gcmeanR(1, nd$day, b[[1]], b[[2]], b[[3]], 0, b[[4]]), tolerance=1e-10)
})

test_that("the asymmetry of each model can be estimated, with the log link", {
  set.seed(4)
  d <- data.frame(day=rep(seq(2, 60, by=4), 3))
  for (m in c("gompertz", "hill", "weibull")){
    mod <- match(m, c("logistic", "gompertz", "hill", "weibull"))
    d$y <- gcTMB:::gcmeanR(mod, d$day, 50, 20, 4, 0, 2) + rnorm(nrow(d), 0, 0.5)
    f <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m, fix=c(b5=FALSE), link=c(b5="log"))
    expect_equal(f$opt$convergence, 0, label=m)
    expect_equal(exp(bcoef(f)[[4]]), 2, tolerance=0.3, label=m)
  }
})
