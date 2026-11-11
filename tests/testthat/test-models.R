# The growth models: mean functions, agreement with nls, and starting values from the linearisation.

test_that("the mean functions are the logistic and gompertz curves", {
  t <- c(-5, 0, 3, 10, 12.5, 30)
  b <- truepars
  expect_equal(gcmean("logistic", t, b), b[1]/(1 + exp(-(t - b[2])/b[3])))
  expect_equal(gcmean("gompertz", t, b), b[1]*exp(-exp(-(t - b[2])/b[3])))
  # b2 is the time of the maximum growth rate: the logistic curve is at half of its asymptote and the
  # gompertz curve at exp(-1) of it
  expect_equal(gcmean("logistic", b[2], b), b[1]/2)
  expect_equal(gcmean("gompertz", b[2], b), b[1]*exp(-1))
  for (m in models){
    slope <- function(t) numDeriv::grad(function(s) gcmean(m, s, b), t)
    tmax <- optimize(slope, c(0, 24), maximum=TRUE, tol=1e-10)$maximum
    expect_equal(tmax, b[2], tolerance=1e-4, label=m)
  }
})

test_that("gaussian fits agree with nls, with the same standard errors (faba bean heights)", {
  for (m in models){
    f <- gcTMB(height ~ DAS, data=FabaHeights25, model=m)
    n <- fabanls(m)
    # to the convergence tolerance of nls
    expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5, label=m)
    # small-sample inference: the covariance of nls
    # nls uses the Gauss-Newton approximation of the Hessian, gcTMB the exact one
    expect_equal(unname(sqrt(diag(vcov(f)))[1:3]), unname(sqrt(diag(vcov(n)))), tolerance=0.05, label=m)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(n)), tolerance=1e-8, label=m)
    expect_equal(residpar(f), summary(n)$sigma, tolerance=1e-6, label=m)
  }
})

test_that("the reserved parameters b4 and b5 are fixed at 0 and not used", {
  d <- simgauss("logistic")
  f <- gcTMB(y ~ day, data=d)
  expect_equal(f$coefnames, paste("b", 1:5, ":(Intercept)", sep=""))
  expect_equal(f$fix, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(f$start[4:5], c(0, 0))
  expect_length(bcoef(f), 3)
  # a value for b4 does not change the curve
  f4 <- gcTMB(y ~ day, data=d, start=c(b4=7, b5=-2))
  expect_equal(f4$opt$objective, f$opt$objective)
  expect_error(gcTMB(y ~ day, rform=b4 ~ 1 | day, data=d), "does not use b4")
  out <- capture.output(print(f))
  expect_false(any(grepl("^Fixed:", out)))
})

test_that("the models and the formula are checked", {
  d <- simgauss("logistic")
  expect_error(gcTMB(y ~ day, data=d, model="loglogistic"), "Unknown model")
  expect_error(gcTMB(y ~ day, data=d, family="gausian"), "Unknown family")
  expect_error(gcTMBstart(y ~ day, data=d, family="gausian"), "Unknown family")
  d$g <- factor(rep(c("a", "b"), length.out=nrow(d)))
  expect_error(gcTMB(y ~ g, data=d), "numeric time")
})

test_that("starting values from the linearised curve are close to the estimates", {
  for (m in models){
    d <- simgauss(m)
    sv <- gcTMBstart(y ~ day, data=d, model=m)
    expect_equal(names(sv), paste("b", 1:5, ":(Intercept)", sep=""))
    f <- gcTMB(y ~ day, data=d, model=m)
    expect_equal(unname(sv[1:3]), unname(bcoef(f)), tolerance=0.1, label=m)
    expect_equal(unname(sv[4:5]), c(0, 0))
  }
})

test_that("the asymptote is profiled when the data do not reach it", {
  # observations only up to the time of the maximum growth rate and a little beyond
  for (m in models){
    d <- simgauss(m, sd=0.5, times=seq(0, 13, 1))
    sv <- gcTMBstart(y ~ day, data=d, model=m)
    expect_gt(sv[[1]], max(d$y))
    expect_equal(unname(sv[1:3]), truepars, tolerance=0.2, label=m)
    f <- gcTMB(y ~ day, data=d, model=m)
    expect_equal(unname(bcoef(f)), unname(coef(nlsfit(m, d))), tolerance=1e-5, label=m)
  }
})

test_that("given starting values are kept and the others are linearised around them", {
  d <- simgauss("logistic", sd=0.5)
  # with b1 given, the regression of qlogis(y/b1) on time gives b2 and b3
  sv <- gcTMBstart(y ~ day, data=d, start=c(b1=50))
  expect_equal(sv[[1]], 50)
  yd <- tapply(d$y, d$day, mean)
  td <- as.numeric(names(yd))
  q <- yd/50
  ok <- q > 0.01 & q < 0.99
  cf <- coef(lm(qlogis(q[ok]) ~ td[ok]))
  expect_equal(unname(sv[2:3]), unname(c(-cf[1]/cf[2], 1/cf[2])), tolerance=1e-8)
  # with b2 also given, b3 from the regression through the origin
  sv <- gcTMBstart(y ~ day, data=d, start=c(b1=50, b2=10))
  expect_equal(sv[[2]], 10)
  z <- qlogis(q[ok]); tt <- td[ok] - 10
  expect_equal(sv[[3]], 1/(sum(tt*z)/sum(tt^2)), tolerance=1e-8)
  # fixed parameters
  f <- gcTMB(y ~ day, data=d, start=c(b1=50), fix=c(b1=TRUE))
  expect_equal(f$start[1], 50)
  expect_length(bcoef(f), 2)
  expect_error(gcTMB(y ~ day, data=d, fix=c(b1=TRUE)), "give the values")
})

test_that("decreasing curves have a negative scale", {
  d <- simgauss("logistic", b=c(50, 10, -2.5))
  f <- gcTMB(y ~ day, data=d)
  expect_equal(unname(bcoef(f)), unname(coef(nlsfit("logistic", d, start=c(50, 10, -2.5)))), tolerance=1e-6)
  expect_lt(bcoef(f)[[3]], 0)
  expect_equal(summary(f)$direction, "decreasing")
  expect_equal(summary(gcTMB(y ~ day, data=simgauss("logistic")))$direction, "increasing")
})

test_that("grid search starting values", {
  d <- simgauss("gompertz")
  sv <- gcTMBstart(y ~ day, data=d, model="gompertz", method="grid", gsn=c(1, 11, 11, 1, 1))
  f <- gcTMB(y ~ day, data=d, model="gompertz")
  # b1 stays at the linearised start; b2 and b3 on the grid of the default ranges
  sl <- gcTMBstart(y ~ day, data=d, model="gompertz")
  expect_equal(sv[[1]], sl[[1]])
  expect_true(sv[[2]] >= -6 && sv[[2]] <= 30)
  expect_true(sv[[3]] >= sl[[3]]/4 && sv[[3]] <= sl[[3]]*4)
  fg <- gcTMB(y ~ day, data=d, model="gompertz", gsn=10)
  expect_equal(coef(fg), coef(f), tolerance=1e-5)
})

test_that("covariates on the parameters agree with nls (faba bean heights)", {
  d <- transform(FabaHeights25, trt=factor(Treatment))
  f <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=d)
  expect_equal(f$coefnames[1:6], c("b1:(Intercept)", "b1:trt0.3", "b1:trt0.5", "b1:trt0.6", "b1:trt0.8", "b2:(Intercept)"))
  # the curve of each treatment, with a common scale: nls with indexed parameters
  fc <- gcTMB(height ~ DAS, fform=b1 + b2 ~ 0 + trt, data=d)
  expect_equal(fc$opt$objective, f$opt$objective, tolerance=1e-8)
  n <- nls(height ~ a[trt]/(1 + exp(-(DAS - c[trt])/s)), data=d, start=list(a=rep(110, 5), c=rep(32, 5), s=11))
  expect_equal(unname(bcoef(fc)), unname(coef(n)), tolerance=1e-5)
  expect_equal(as.numeric(logLik(fc)), as.numeric(logLik(n)), tolerance=1e-8)
})

test_that("poisson fits maximise the poisson likelihood (faba bean pods)", {
  d <- FabaCounts23
  for (m in models){
    f <- gcTMB(pods ~ DAS, data=d, model=m, family="poisson")
    b <- unname(bcoef(f))
    nll <- function(b){
      mu <- gcmean(m, d$DAS, b)
      if (any(mu < 0)) Inf else -sum(dpois(d$pods, mu, log=TRUE))
    }
    # the log-likelihood of the counts at the estimates, which are a maximum
    expect_equal(as.numeric(logLik(f)), -nll(b), tolerance=1e-10, label=m)
    expect_equal(numDeriv::grad(nll, b), c(0, 0, 0), tolerance=1e-4, label=m)
    # no better optimum from another start
    o <- optim(b*c(1.1, 0.95, 1.2), nll, control=list(reltol=1e-14, maxit=5000))
    expect_gte(o$value, f$opt$objective - 1e-6)
    # the covariance is the inverse of the Hessian of the negative log-likelihood
    # (numerical derivatives with steps of 1% of the parameters; the default 10% are too coarse)
    expect_equal(unname(vcov(f)), solve(numDeriv::hessian(nll, b, method.args=list(d=0.01))), tolerance=1e-6, label=m)
    # no residual variance parameter, and no small-sample inference
    expect_length(coef(f), 3)
    expect_null(df.residual(f))
  }
})

test_that("poisson models with covariates and random effects (faba bean pods)", {
  d <- FabaCounts23
  f0 <- gcTMB(pods ~ DAS, data=d, family="poisson")
  f1 <- gcTMB(pods ~ DAS, fform=b1 ~ Genotype, data=d, family="poisson")
  f2 <- gcTMB(pods ~ DAS, fform=b1 ~ Genotype, rform=b1 ~ 1 | Tube_ID, data=d, family="poisson")
  expect_length(bcoef(f1), 8)
  expect_gt(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_gt(as.numeric(logLik(f2)), as.numeric(logLik(f1)))
  # the plants vary in their final number of pods
  expect_gt(gcTMB:::recov(f2)[[1]]$sd, 1)
  expect_true(f2$sdrl$pdHess)
  s <- summary(f2)
  expect_equal(s$ngroups, c(Tube_ID=22))
  # Pearson dispersion close to 1 conditional on the random effects
  expect_lt(abs(s$dispersion[["ratio"]] - 1), 0.3)
  # the asymptote of each genotype from predict
  nd <- data.frame(Genotype=sort(unique(d$Genotype)))
  pp <- predict(f1, newdata=nd, type="parameters")
  expect_equal(unname(pp[, "b1"]), unname(bcoef(f1)[1] + c(0, bcoef(f1)[2:6])))
})
