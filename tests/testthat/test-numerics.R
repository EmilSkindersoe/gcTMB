# Internal numerical checks: automatic derivatives against finite differences, and invariances of the
# fitted models.

skip_if_not_installed("numDeriv")

# a point away from the optimum, where the gradient is not zero
offopt <- function(fit) fit$opt$par * 1.05 + 0.01

test_that("gradients and hessians agree with finite differences for the models and links", {
  for (m in models){
    for (lk in list(NULL, c(b3="inverse"))){
      fit <- gcTMB(y ~ day, data=simgauss(m), model=m, link=lk)
      p <- offopt(fit)
      expect_equal(as.vector(fit$obj$gr(p)), numDeriv::grad(fit$obj$fn, p), tolerance=1e-6, label=m)
      expect_equal(unname(fit$obj$he(p)), numDeriv::jacobian(fit$obj$gr, p), tolerance=1e-6, label=m)
    }
  }
})

test_that("gradients agree with finite differences for the families", {
  set.seed(5)
  d <- data.frame(day=rep(days, each=4), n=20)
  mu <- gcmean("logistic", d$day, c(0.9, 10, 2.5))
  d$yb <- rbinom(nrow(d), d$n, mu)
  d$yp <- rpois(nrow(d), 40*mu)
  d$yn <- rnbinom(nrow(d), size=5, mu=40*mu)
  d$yr <- rbeta(nrow(d), 30*mu, 30*(1 - mu))
  d$yl <- 40*mu*exp(rnorm(nrow(d), 0, 0.1)) + 1
  fits <- list(
    binomial=gcTMB(cbind(yb, n - yb) ~ day, data=d, family="binomial"),
    poisson=gcTMB(pods ~ DAS, data=FabaCounts23, family="poisson", model="gompertz"),
    nbinom2=gcTMB(yn ~ day, data=d, family="nbinom2"),
    beta=gcTMB(yr ~ day, data=d, family="beta"),
    lognormal=gcTMB(yl ~ day, data=d, family="lognormal")
  )
  for (nm in names(fits)){
    fit <- fits[[nm]]
    p <- fit$opt$par * 1.01
    expect_equal(as.vector(fit$obj$gr(p)), numDeriv::grad(fit$obj$fn, p), tolerance=1e-5, label=nm)
  }
})

test_that("gradients of the Laplace approximation agree with finite differences", {
  d <- simre()
  fits <- list(
    single=gcTMB(y ~ day, rform=b2 ~ 1 | id, data=d),
    correlated=gcTMB(y ~ day, rform=b1 + b2 ~ 1 | id, data=d),
    independent=gcTMB(y ~ day, rform=list(b1 ~ 1 | id, b2 ~ 1 | id), data=d),
    rate=gcTMB(y ~ day, rform=b3 ~ 1 | id, data=d, link=c(b3="inverse"))
  )
  for (nm in names(fits)){
    fit <- fits[[nm]]
    p <- offopt(fit)
    expect_equal(as.vector(fit$obj$gr(p)), numDeriv::grad(fit$obj$fn, p), tolerance=1e-5, label=nm)
  }
})

test_that("weights are equivalent to replicated observations", {
  d <- simgauss("logistic")
  w <- rep(1:2, length.out=nrow(d))
  ml <- list(smallsample=FALSE)
  f1 <- gcTMB(y ~ day, data=d, weights=w, control=ml)
  f2 <- gcTMB(y ~ day, data=d[rep(seq_len(nrow(d)), w), ], control=ml)
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-8)
  expect_equal(coef(f1), coef(f2), tolerance=1e-5)
  expect_equal(vcov(f1), vcov(f2), tolerance=1e-4)
})

test_that("the order of the observations does not matter", {
  d <- simre()
  set.seed(8)
  dp <- d[sample(nrow(d)), ]
  f1 <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=d)
  f2 <- gcTMB(y ~ day, rform=b1 ~ 1 | id, data=dp)
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-8)
  expect_equal(coef(f1), coef(f2), tolerance=1e-5)
})

test_that("contrasts of the fixed effects do not change the fit", {
  d <- rbind(cbind(simgauss("logistic", seed=1), trt="a"),
             cbind(simgauss("logistic", seed=2, b=c(55, 11, 2.5)), trt="b"),
             cbind(simgauss("logistic", seed=3, b=c(45, 9, 2.5)), trt="c"))
  d$trt <- factor(d$trt)
  f1 <- gcTMB(y ~ day, fform=b1 + b2 ~ trt, data=d)
  d2 <- d
  contrasts(d2$trt) <- contr.sum(3)
  f2 <- gcTMB(y ~ day, fform=b1 + b2 ~ trt, data=d2)
  f3 <- gcTMB(y ~ day, fform=b1 + b2 ~ 0 + trt, data=d)
  expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-7)
  expect_equal(f1$opt$objective, f3$opt$objective, tolerance=1e-7)
  expect_equal(predict(f1), predict(f3), tolerance=1e-5)
})

test_that("shifting and rescaling time changes b2 and b3 only", {
  for (m in models){
    d <- simgauss(m)
    f1 <- gcTMB(y ~ day, data=d, model=m)
    d2 <- transform(d, day=day*24 + 100)
    f2 <- gcTMB(y ~ day, data=d2, model=m)
    expect_equal(f1$opt$objective, f2$opt$objective, tolerance=1e-7)
    b <- bcoef(f1)
    expect_equal(unname(bcoef(f2)), unname(c(b[1], b[2]*24 + 100, b[3]*24)), tolerance=1e-5, label=m)
    expect_equal(unname(coef(timeto(f2, 90))), unname(coef(timeto(f1, 90))*24 + 100), tolerance=1e-5)
  }
})

test_that("the residual parameter starts at its maximum likelihood value for the starting curve", {
  d <- simgauss("logistic")
  w <- rep(1:3, length.out=nrow(d))
  f <- gcTMB(y ~ day, data=d, weights=w)
  s <- f$start
  mu0 <- gcmean("logistic", d$day, s[1:3])
  expect_equal(f$obj$par[["log_sigma"]], log(sqrt(sum(w*(d$y - mu0)^2)/sum(w))), tolerance=1e-6)
  # on the scale of the data
  f4 <- gcTMB(y ~ day, data=transform(d, y=y*1e4), weights=w)
  expect_equal(f4$obj$par[["log_sigma"]] - f$obj$par[["log_sigma"]], log(1e4), tolerance=1e-6)
})
