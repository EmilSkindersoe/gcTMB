# The additive parameter b6: the curve b6 + f(t), fixed at b6 = 0 by default

test_that("b6 is added to the curves and fixed at 0 by default", {
  t <- c(0, 5, 10, 20)
  for (mod in c(1:6, 8:9)){
    expect_equal(gcTMB:::gcmeanR(mod, t, 50, 10, 3, 0, 1, f6=1.3), 1.3 + gcTMB:::gcmeanR(mod, t, 50, 10, 3, 0, 1), label=mod)
  }
  # the linear model has no b2-b5: b6 is its intercept
  expect_equal(gcTMB:::unusedpars(7), 2:5)
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  expect_equal(f$coefnames[6], "b6:(Intercept)")
  expect_true(f$fix[6])
  expect_equal(f$start[6], 0)
  expect_false(any(grepl("b6", names(coef(f)))))
  expect_equal(colnames(predict(f, type="parameters")), c("b1", "b2", "b3"))
})

test_that("a fixed and an estimated b6 shift the curve, and timeto and maxrate do not change", {
  d <- simgauss("logistic")
  f <- gcTMB(y ~ day, data=d)
  d20 <- transform(d, y=y + 20)
  # b6 fixed at the known shift: the same fit
  f20 <- gcTMB(y ~ day, data=d20, start=c(b6=20))
  expect_true(f20$fix[6])
  expect_equal(bcoef(f20)[1:3], bcoef(f)[1:3], tolerance=1e-6)
  expect_equal(predict(f20), predict(f) + 20, tolerance=1e-6)
  expect_equal(coef(timeto(f20, level=c(10, 50, 90))), coef(timeto(f, level=c(10, 50, 90))), tolerance=1e-6)
  expect_equal(vcov(timeto(f20, level=50)), vcov(timeto(f, level=50)), tolerance=1e-4)
  expect_equal(coef(maxrate(f20)), coef(maxrate(f)), tolerance=1e-6)
  # absolute times include b6
  expect_equal(coef(timeto(f20, level=40, type="absolute")), coef(timeto(f, level=20, type="absolute")),
               tolerance=1e-6, ignore_attr=TRUE)
  # a fixed b6 is listed with the fixed coefficients
  expect_true(any(grepl("b6:\\(Intercept\\) = 20", capture.output(print(f20)))))
  expect_equal(colnames(predict(f20, type="parameters")), c("b1", "b2", "b3", "b6"))

  # b6 estimated: with observations before the growth (b4 = 0 truncates), b6 is identified
  set.seed(2)
  d <- data.frame(day=rep(seq(-10, 24, 2), each=3))
  d$y <- 20 + gcTMB:::gcmeanR(1, d$day, 50, 10, 2.5, 0, trunc=TRUE) + rnorm(nrow(d), 0, 1.5)
  fe <- gcTMB(y ~ day, data=d, fix=c(b6=FALSE))
  expect_false(fe$fix[6])
  expect_equal(fe$opt$convergence, 0)
  expect_equal(coef(fe)[["b6:(Intercept)"]], 20, tolerance=0.05)
  n <- nls(y ~ b6 + ifelse(day > 0, b1*(plogis((day - b2)/b3) - plogis(-b2/b3))/(1 - plogis(-b2/b3)), 0), data=d,
           start=list(b1=50, b2=10, b3=2.5, b6=20))
  expect_equal(unname(coef(fe)[grepl("^b", names(coef(fe)))]), unname(coef(n)), tolerance=1e-5)
  expect_equal(sqrt(vcov(fe)[4, 4]), sqrt(vcov(n)[4, 4]), tolerance=2e-2)
  # the log link keeps b6 positive
  fl <- gcTMB(y ~ day, data=d, fix=c(b6=FALSE), link=c(b6="log"))
  expect_equal(exp(coef(fl)[["b6:(Intercept)"]]), coef(fe)[["b6:(Intercept)"]], tolerance=1e-5)
})

test_that("b6 for a height-diameter model starting from breast height", {
  set.seed(4)
  d <- data.frame(dbh=runif(150, 5, 60))
  d$h <- 1.3 + gcTMB:::gcmeanR(8, d$dbh, 30, 8, 12, 0) + rnorm(150, 0, 1)
  f <- gcTMB(h ~ dbh, data=d, model="chapmanrichards", start=c(b6=1.3))
  n <- nls(h ~ 1.3 + b1*(1 - exp(-dbh/b3))^exp(b2/b3), data=d, start=list(b1=30, b2=8, b3=12))
  expect_equal(unname(bcoef(f)[1:3]), unname(coef(n)), tolerance=1e-5)
  p <- plot(f, which=2, draw=FALSE)
  expect_equal(p$time$labels$title, "Residuals vs dbh")
})

test_that("covariates and random effects on b6", {
  set.seed(5)
  d <- simgroups(ng=10)
  u <- rnorm(10, 0, 3)
  d$y <- 10 + u[d$id] + gcTMB:::gcmeanR(1, d$day, 50, 10, 2.5, 0, trunc=TRUE) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, rform=b6 ~ 1 | id, data=d, start=c(b6=10), fix=c(b6=FALSE))
  expect_equal(f$opt$convergence, 0)
  expect_equal(unname(gcTMB:::recov(f)[[1]]$sd), 3, tolerance=0.5)
  d$g <- factor(as.integer(d$id) %% 2)
  fg <- gcTMB(y ~ day, fform=b6 ~ g, data=d, start=c(b6=10), fix=c(b6=FALSE))
  expect_equal(fg$coefnames[6:7], c("b6:(Intercept)", "b6:g1"))
})

test_that("start, fix, bounds and curves without b6 are still accepted", {
  d <- simgauss("logistic")
  f <- gcTMB(y ~ day, data=d)
  # unnamed vectors for b1-b5
  f5 <- gcTMB(y ~ day, data=d, start=c(50, 10, 2.5, 0, 1), fix=c(FALSE, FALSE, FALSE, TRUE, TRUE),
              lower=c(0, -Inf, 0, -Inf, -Inf), upper=rep(Inf, 5))
  expect_equal(f5$fix, f$fix)
  expect_equal(bcoef(f5), bcoef(f), tolerance=1e-6)
  expect_equal(length(gcTMBstart(y ~ day, data=d, method="grid", gsn=c(1, 5, 5, 1, 1))), 6)
  # a list of five right-hand side formulas
  ff <- gcTMB(y ~ day, fform=list(~ 1, ~ 1, ~ 1, ~ 1, ~ 1), data=d)
  expect_equal(bcoef(ff), bcoef(f), tolerance=1e-6)
  # curves with 5 rows and columns
  expect_equal(coef(timeto(f, curves=diag(5))), coef(timeto(f)))
  expect_equal(coef(maxrate(f, curves=diag(5))), coef(maxrate(f)))
  # pairwise accepts b6
  expect_error(pairwise(f, newdata=data.frame(x=1:2), parm="b7"), "b1, b2, b3, b4, b5 and b6")
})
