# Times at which the curves reach a level of their asymptote.

test_that("times are b2 + b3 Finv(p), with standard errors of the delta method", {
  for (m in models){
    f <- gcTMB(y ~ day, data=simgauss(m), model=m)
    b <- bcoef(f)
    p <- c(10, 50, 90)
    tt <- timeto(f, respLev=p)
    expect_s3_class(tt, "timetoTMB")
    expect_equal(names(coef(tt)), paste(p, ":lf", 1:3, sep=""))
    expect_equal(unname(coef(tt)), b[[2]] + b[[3]]*Finv[[m]](p/100), tolerance=1e-10, label=m)
    V <- vcov(f)[1:3, 1:3]
    se <- sapply(p, function(pp){
      g <- c(0, 1, Finv[[m]](pp/100))
      sqrt(drop(t(g) %*% V %*% g))
    })
    expect_equal(unname(sqrt(diag(vcov(tt)))), se, tolerance=1e-8, label=m)
    # the curve reaches the level at the time
    expect_equal(predict(f, newdata=data.frame(day=coef(tt))), b[[1]]*p/100, tolerance=1e-8)
  }
  # the logistic curve reaches half of its asymptote at b2
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  expect_equal(unname(coef(timeto(f))), bcoef(f)[[2]])
  expect_equal(unname(vcov(timeto(f))[1, 1]), vcov(f)[2, 2])
})

test_that("times do not depend on the link of b3", {
  d <- simgauss("gompertz")
  f <- gcTMB(y ~ day, data=d, model="gompertz")
  fi <- gcTMB(y ~ day, data=d, model="gompertz", link=c(b3="inverse"))
  expect_equal(coef(timeto(fi, c(5, 95))), coef(timeto(f, c(5, 95))), tolerance=1e-6)
  expect_equal(vcov(timeto(fi, c(5, 95))), vcov(timeto(f, c(5, 95))), tolerance=1e-4)
  expect_equal(suppressMessages(confint(timeto(fi, 95))), suppressMessages(confint(timeto(f, 95))), tolerance=1e-4)
})

test_that("absolute times: times at which the curve reaches a response", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  b <- bcoef(f)
  tt <- timeto(f, respLev=c(10, 40), type="absolute")
  expect_equal(names(coef(tt)), c("y=10:lf1", "y=40:lf2"))
  expect_equal(unname(coef(tt)), b[[2]] + b[[3]]*qlogis(c(10, 40)/b[[1]]), tolerance=1e-10)
  # the standard errors include the asymptote
  g <- function(bb) bb[2] + bb[3]*qlogis(40/bb[1])
  gr <- numDeriv::grad(g, unname(b))
  expect_equal(sqrt(vcov(tt)[2, 2]), sqrt(drop(t(gr) %*% vcov(f)[1:3, 1:3] %*% gr)), tolerance=1e-6)
  # responses beyond the asymptote are not reached
  expect_warning(tt <- timeto(f, respLev=c(20, 80), type="absolute"), "do not reach")
  expect_true(is.na(coef(tt)[2]))
  expect_error(timeto(f, respLev=100), "between 0 and 100")
})

test_that("timeto with newdata gives the times of the curve of each row", {
  d <- rbind(cbind(simgauss("logistic", seed=1), trt="a"),
             cbind(simgauss("logistic", seed=2, b=c(60, 12, 3)), trt="b"))
  d$trt <- factor(d$trt)
  f <- gcTMB(y ~ day, fform=list(b1 ~ trt, b2 ~ trt, b3 ~ trt), data=d)
  tt <- timeto(f, c(50, 90), newdata=data.frame(trt=c("a", "b")))
  expect_equal(names(coef(tt)), c("50:trt=a", "50:trt=b", "90:trt=a", "90:trt=b"))
  b <- bcoef(f)
  expect_equal(unname(coef(tt)[1:2]), unname(c(b[3], b[3] + b[4])))
  # the same as linfct
  L <- matrix(0, 5, 8)
  L[1, 1] <- L[2, 3] <- L[3, 5] <- 1
  L[1, 2] <- L[2, 4] <- L[3, 6] <- 1
  expect_equal(unname(coef(timeto(f, 90, linfct=L))), unname(coef(tt)[4]))
  expect_error(timeto(f), "8 model coefficients")
  expect_error(timeto(f, newdata=data.frame(x=1)), "covariates of fform")
  expect_error(timeto(f, linfct=L, newdata=data.frame(trt="a")), "either linfct or newdata")
})

test_that("times of fitted groups and of the population", {
  d <- simre()
  f <- gcTMB(y ~ day, rform=b1 + b2 ~ 1 | id, data=d)
  # population: the fixed effects
  tp <- timeto(f, 50, random=FALSE)
  expect_equal(unname(coef(tp)), bcoef(f)[[2]])
  expect_equal(unname(vcov(tp)[1, 1]), vcov(f)[2, 2])
  # a new group: the variance of the random effects of b2 is added
  tn <- timeto(f, 50)
  vc <- gcTMB:::recov(f)[[1]]
  expect_equal(unname(vcov(tn)[1, 1]), vcov(f)[2, 2] + vc$cov[2, 2], tolerance=1e-8)
  # fitted groups: the curves with their predicted random effects
  tg <- timeto(f, 50, newdata=data.frame(id=c("1", "2")))
  pp <- predict(f, newdata=data.frame(id=c("1", "2")), type="parameters")
  expect_equal(unname(coef(tg)), unname(pp[, "b2"]), tolerance=1e-8)
  expect_true(all(tg$conditional))
  # profile intervals are not available for fitted groups
  expect_message(ci <- confint(tg), "Wald intervals")
  expect_equal(dim(ci), c(2, 2))
})

test_that("residual variance only for the gaussian family", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  t0 <- timeto(f, 90)
  t1 <- timeto(f, 90, residual=TRUE)
  expect_equal(coef(t1), coef(t0))
  expect_gt(vcov(t1)[1, 1], vcov(t0)[1, 1])
  set.seed(3)
  dp <- data.frame(day=rep(days, each=3))
  dp$y <- rpois(nrow(dp), gcmean("logistic", dp$day, truepars))
  fp <- gcTMB(y ~ day, data=dp, family="poisson")
  expect_error(timeto(fp, residual=TRUE), "gaussian")
})

test_that("print returns the times invisibly", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  tt <- timeto(f)
  out <- capture.output(res <- withVisible(print(tt)))
  expect_false(res$visible)
  expect_identical(res$value, tt)
  expect_true(any(grepl("^Times to reach the percentages of the asymptote", out)))
})
