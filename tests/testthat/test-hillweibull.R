# The hill and weibull models on the log time axis: f(t) = b1 / (1 + (b2/t)^b3) and
# f(t) = b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3)), zero up to t = 0.

lt <- c("hill", "weibull")
ltpars <- list(hill=c(50, 10, 4), weibull=c(50, 10, 3))

test_that("the curves are the hill and weibull functions, zero up to t = 0", {
  t <- c(-1, 0, 0.5, 5, 10, 20, 50)
  b <- c(50, 10, 3)
  tp <- pmax(t, 0)
  expect_equal(gcTMB:::gcmeanR(3, t, b[1], b[2], b[3], 0, 1), ifelse(t > 0, b[1]/(1 + (b[2]/tp)^b[3]), 0))
  expect_equal(gcTMB:::gcmeanR(4, t, b[1], b[2], b[3], 0, 1), ifelse(t > 0, b[1]*(1 - exp(-(b[3] - 1)/b[3]*(tp/b[2])^b[3])), 0))
  # b2 is the time of half the asymptote (hill) and the inflection point (weibull)
  expect_equal(gcTMB:::gcmeanR(3, 10, 1, 10, 3, 0, 1), 0.5)
  w <- function(t) gcTMB:::gcmeanR(4, t, 50, 10, 3, 0, 1)
  expect_equal(optimize(function(t) -(w(t + 1e-4) - w(t - 1e-4)), c(1, 30), tol=1e-8)$minimum, 10, tolerance=1e-4)
  # times to reach a level are the inverse of the curves, also truncated
  q <- c(0.1, 0.5, 0.9)
  for (mod in 3:4){
    for (trunc in c(FALSE, TRUE)){
      tt <- gcTMB:::gctimeR(mod, q, 10, 3, 4, trunc=trunc)
      expect_equal(gcTMB:::gcmeanR(mod, tt, 1, 10, 3, 4, 1, trunc=trunc), q, tolerance=1e-12)
    }
    # a threshold up to zero does not truncate the curve
    expect_equal(gcTMB:::gcmeanR(mod, t, 50, 10, 3, -2, 1, trunc=TRUE), gcTMB:::gcmeanR(mod, t, 50, 10, 3, 0, 1))
  }
  expect_equal(gcTMB:::gctimeR(3, 0.5, 10, 3), 10)
})

# the weibull shape b3 on its own scale, instead of its default link log(b3 - 1)
lkid <- function(m) if (m == "weibull") c(b3="identity")

test_that("hill and weibull fits agree with nls", {
  for (m in lt){
    d <- simgauss(m, b=ltpars[[m]])
    f <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m, link=lkid(m))
    expect_equal(f$opt$convergence, 0, label=m)
    n <- nlsfit(m, d, start=ltpars[[m]])
    expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5, label=m)
    # observed Hessian (gcTMB) and its Gauss-Newton approximation (nls)
    expect_equal(unname(sqrt(diag(vcov(f)))[1:3]), unname(sqrt(diag(vcov(n)))), tolerance=0.05, label=m)
    expect_equal(as.vector(predict(f)), as.vector(fitted(n)), tolerance=1e-6, label=m)
    expect_match(grep("^Model:", capture.output(print(f)), value=TRUE), m)
  }
})

test_that("starting values of the hill and weibull models", {
  for (m in lt){
    d <- simgauss(m, b=ltpars[[m]])
    sv <- gcTMBstart(y ~ day, truncated=FALSE, data=d, model=m, link=lkid(m))
    expect_true(all(is.finite(sv[1:3])))
    expect_equal(unname(sv[1:3]), ltpars[[m]], tolerance=0.3, label=m)
    # with b2 given, b3 from the linearised curve
    sv2 <- gcTMBstart(y ~ day, truncated=FALSE, data=d, model=m, start=c(b2=10), link=lkid(m))
    expect_equal(unname(sv2[3]), ltpars[[m]][3], tolerance=0.3, label=m)
  }
  # the weibull shape starts above 1
  expect_gt(gcTMBstart(height ~ DAS, truncated=FALSE, data=FabaHeights25, model="weibull", link=c(b3="identity"))[[3]], 1)
})

test_that("times of the hill and weibull curves, with standard errors and profiles", {
  for (m in lt){
    d <- simgauss(m, b=ltpars[[m]])
    f <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m, link=lkid(m))
    b <- bcoef(f)
    p <- c(10, 50, 90)
    tt <- timeto(f, level=p)
    expect_equal(unname(coef(tt)), gcTMB:::gctimeR(match(m, allmodels), p/100, b[[2]], b[[3]]), tolerance=1e-10)
    expect_equal(unname(predict(f, newdata=data.frame(day=coef(tt)))), b[[1]]*p/100, tolerance=1e-8)
    # delta method, by numerical derivatives
    g <- numDeriv::jacobian(function(x) gcTMB:::gctimeR(match(m, allmodels), p/100, x[2], x[3]), b)
    expect_equal(unname(vcov(tt)), unname(g %*% vcov(f)[1:3, 1:3] %*% t(g)), tolerance=1e-5, label=m)
    ci <- confint(tt)
    expect_true(all(is.finite(ci)))
    expect_true(all(ci[, 1] < coef(tt) & coef(tt) < ci[, 2]))
  }
  # hill: t50 is b2
  f <- gcTMB(y ~ day, truncated=FALSE, data=simgauss("hill", b=ltpars$hill), model="hill")
  expect_equal(unname(coef(timeto(f, level=50))), bcoef(f)[[2]], tolerance=1e-10)
})

test_that("truncated hill and weibull curves recover the threshold", {
  for (m in lt){
    set.seed(5)
    d <- data.frame(day=rep(0:30, each=4))
    d$y <- gcTMB:::gcmeanR(match(m, allmodels), d$day, 50, 12, 3, 6, 1, trunc=TRUE) + rnorm(nrow(d), 0, 1)
    f <- gcTMB(y ~ day, data=d, model=m, fix=c(b4=FALSE))
    expect_equal(f$opt$convergence, 0, label=m)
    expect_equal(coef(f)[["b4:(Intercept)"]], 6, tolerance=0.1, label=m)
    expect_lt(AIC(f), AIC(gcTMB(y ~ day, truncated=FALSE, data=d, model=m)))
    expect_equal(as.vector(f$obj$gr(f$opt$par)), numDeriv::grad(f$obj$fn, f$opt$par), tolerance=1e-5)
    tt <- timeto(f, level=c(10, 50))
    expect_equal(unname(predict(f, newdata=data.frame(day=coef(tt)))), coef(f)[["b1:(Intercept)"]]*c(0.1, 0.5), tolerance=1e-8)
    # profiles of the times: b2 found numerically on the log scale
    ci <- confint(tt)
    expect_true(all(is.finite(ci)))
    expect_true(all(ci[, 1] < coef(tt) & coef(tt) < ci[, 2]))
  }
})

test_that("hill and weibull models with covariates, random effects, links and other families", {
  FabaHeights25$trt <- factor(FabaHeights25$Treatment)
  for (m in lt){
    f <- gcTMB(height ~ DAS, truncated=FALSE, data=FabaHeights25, model=m)
    expect_equal(f$opt$convergence, 0, label=m)
    fl <- gcTMB(height ~ DAS, truncated=FALSE, data=FabaHeights25, model=m, link=c(b2="log", b3="log"))
    expect_equal(fl$opt$objective, f$opt$objective, tolerance=1e-7, label=m)
    fr <- gcTMB(height ~ DAS, truncated=FALSE, fform=b1 ~ trt, rform=b1 ~ 1 | plant_ID, data=FabaHeights25, model=m)
    expect_equal(fr$opt$convergence, 0, label=m)
    fp <- gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson", model=m, fix=c(b4=FALSE))
    expect_equal(fp$opt$convergence, 0, label=m)
    expect_true(is.finite(logLik(fp)))
  }
  expect_error(gcTMB(height ~ DAS, truncated=FALSE, data=FabaHeights25, model="hill", link=c(b3="inverse")), "not available for b3")
})
