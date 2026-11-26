# maxrate: the time of the maximum growth rate, the rate and the fraction of the asymptote reached then

# the maximum of the slope of a curve on a fine grid
gridmax <- function(mod, b, trunc=FALSE){
  t <- seq(max(0.01, b[4] - 5), 150, by=0.001)
  f <- gcTMB:::gcmeanR(mod, t, b[1], b[2], b[3], b[4], b[5], trunc=trunc)
  r <- diff(f)/0.001
  i <- which.max(r)
  c(tmax=t[i] + 0.0005, rate=r[i], fraction=mean(f[i + 0:1])/b[1])
}

test_that("tmax, rate and fraction are the maximum of the slope of the curves", {
  pars <- list(c(50, 30, 5), c(50, 30, 5), c(50, 30, 4), c(50, 30, 3))
  for (mod in 1:4){
    for (a in c(0.5, 1, 2.5)){
      b <- c(pars[[mod]], 0, a)
      expect_equal(gcTMB:::maxrateR(mod, b), gridmax(mod, b), tolerance=1e-3, label=paste(mod, a))
    }
  }
  # closed forms of the symmetric logistic and gompertz curves
  expect_equal(gcTMB:::maxrateR(1, c(50, 30, 5, 0, 1)), c(tmax=30, rate=50/20, fraction=0.5))
  expect_equal(gcTMB:::maxrateR(2, c(50, 30, 5, 0, 1)), c(tmax=30, rate=50/(exp(1)*5), fraction=exp(-1)))
})

test_that("truncated curves grow fastest at the inflection point or just after the threshold", {
  for (mod in 1:4){
    # threshold before the inflection point
    b <- c(50, 30, 4, 15, 1.5)
    expect_equal(gcTMB:::maxrateR(mod, b, trunc=TRUE), gridmax(mod, b, trunc=TRUE), tolerance=1e-3, label=mod)
    # threshold after the inflection point: the rate just after the threshold
    b <- c(50, 30, 4, 40, 1.5)
    mr <- gcTMB:::maxrateR(mod, b, trunc=TRUE)
    expect_equal(mr[["tmax"]], 40)
    expect_equal(mr[["fraction"]], 0)
    h <- 1e-6
    slope <- (gcTMB:::gcmeanR(mod, 40 + 2*h, 50, 30, 4, 40, 1.5, trunc=TRUE) - gcTMB:::gcmeanR(mod, 40 + h, 50, 30, 4, 40, 1.5, trunc=TRUE))/h
    expect_equal(mr[["rate"]], slope, tolerance=1e-4, label=mod)
  }
})

test_that("hill and weibull curves without an inflection point give NA", {
  expect_true(all(is.na(gcTMB:::maxrateR(3, c(50, 30, 2, 0, 0.4)))))
  expect_true(all(is.na(gcTMB:::maxrateR(4, c(50, 30, 2, 0, 0.4)))))
})

test_that("maxrate has standard errors of the delta method and Wald, profile and bootstrap intervals", {
  d <- simgauss("logistic")
  f <- gcTMB(y ~ day, truncated=FALSE, data=d, fix=c(b5=FALSE), link=c(b5="log"), start=c(b5=0))
  mr <- maxrate(f)
  expect_equal(names(coef(mr)), c("tmax", "rate", "fraction"))
  b <- bcoef(f)
  val <- function(p) gcTMB:::maxrateR(1, c(p[1:3], 0, exp(p[4])))
  expect_equal(unname(coef(mr)), unname(val(b)))
  J <- numDeriv::jacobian(val, unname(b))
  expect_equal(unname(vcov(mr)), unname(J %*% vcov(f)[1:4, 1:4] %*% t(J)), tolerance=1e-5)
  w <- confint(mr, method="wald")
  expect_equal(unname(w[, 2] - w[, 1]), unname(2*qt(0.975, f$rdf)*sqrt(diag(vcov(mr)))))
  p <- suppressMessages(confint(mr, parm=c("tmax", "rate"), method="profile"))
  expect_equal(p, w[1:2, ], tolerance=0.05)
  # the profile is the likelihood with the quantity fixed: refit with b2 or b1 at a limit
  lim <- p[1, 1]
  fl <- function(b2) gcTMB(y ~ day, truncated=FALSE, data=d, fix=c(b2=TRUE, b5=FALSE), link=c(b5="log"), start=c(b2=b2, b5=0))$opt$objective
  # for the Richards curve, tmax is b2
  expect_equal(2*(fl(lim) - f$opt$objective), gcTMB:::profcutoff(0.95, f), tolerance=1e-3)
  bt <- bootstrap(f, nsim=20, seed=1, progress=FALSE)
  bci <- confint(mr, boot=bt)
  expect_true(all(bci[, 1] <= coef(mr) + 1e-8 & bci[, 2] >= coef(mr) - 1e-8))
})

test_that("maxrate for curves of newdata and for truncated fits", {
  d <- simgauss("logistic")
  d$g <- factor(rep(c("a", "b"), length.out=nrow(d)))
  f <- gcTMB(y ~ day, truncated=FALSE, fform=b1 ~ g, data=d)
  mr <- maxrate(f, newdata=data.frame(g=c("a", "b")))
  b <- bcoef(f)
  expect_equal(unname(coef(mr)[c("rate:g=a", "rate:g=b")]), c(b[[1]], b[[1]] + b[[2]])/(4*b[[4]]))
  ft <- gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", fix=c(b4=FALSE))
  bt <- coef(ft)[1:4]
  expect_equal(unname(coef(maxrate(ft))), unname(gcTMB:::maxrateR(1, c(bt, 1), trunc=TRUE)))
  # a hill curve with b3 b5 < 1, by a linear function that scales b5 = 1 down to 0.01
  fh <- gcTMB(Height ~ DAS, truncated=FALSE, data=FabaHeights25, model="hill")
  expect_warning(mh <- maxrate(fh, curves=diag(c(1, 1, 1, 1, 0.01))), "no inflection point")
  expect_true(all(is.na(coef(mh))))
})

test_that("maxrate gives the curves of fitted groups, with their predicted random effects", {
  d <- FabaHeights25
  # without truncation: tmax = b2, rate = b1/(4 b3), fraction 1/2 of the plant's own curve
  f <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, truncated=FALSE, data=d)
  nd <- data.frame(Plant_ID=c(6, 7))
  mr <- maxrate(f, newdata=nd)
  expect_equal(names(coef(mr))[1:3], c("tmax:Plant_ID=6", "rate:Plant_ID=6", "fraction:Plant_ID=6"))
  p <- predict(f, newdata=nd, type="parameters")
  expect_equal(unname(coef(mr)[c(1, 4)]), unname(p[, "b2"]), tolerance=1e-6)
  expect_equal(unname(coef(mr)[c(2, 5)]), unname(p[, "b1"]/(4*p[, "b3"])), tolerance=1e-6)
  expect_equal(unname(coef(mr)[c(3, 6)]), c(0.5, 0.5))
  # the plants differ from the population curve; random=FALSE gives the population curve
  pop <- coef(maxrate(f))
  expect_false(isTRUE(all.equal(unname(coef(mr)[1:3]), unname(pop))))
  expect_equal(unname(coef(maxrate(f, newdata=nd[1, , drop=FALSE], random=FALSE))), unname(pop))
  # the standard error of tmax is that of the time to 50% of the plant's curve, its inflection point
  tt <- timeto(f, 50, newdata=nd)
  expect_equal(unname(sqrt(diag(vcov(mr)))[c(1, 4)]), unname(sqrt(diag(vcov(tt)))), tolerance=1e-4)
  # profile intervals are not available for fitted groups: Wald intervals, with a message
  expect_message(ci <- confint(mr, method="profile"), "population curves only")
  expect_equal(ci, confint(mr))
  expect_error(profile(mr, parm="rate:Plant_ID=6"), "population curves")
  # a level that is not in the fitted data: the population curve, with a warning
  expect_warning(mu <- maxrate(f, newdata=data.frame(Plant_ID=999)), "not in the fitted data")
  expect_equal(unname(coef(mu)), unname(pop))
})

test_that("maxrate of a fitted group of a truncated logistic curve is rescaled by the truncation", {
  f <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25)
  nd <- data.frame(Plant_ID=6)
  mr <- maxrate(f, newdata=nd)
  p <- predict(f, newdata=nd, type="parameters")
  b1 <- p[, "b1"]; b2 <- p[, "b2"]; b3 <- p[, "b3"]; b4 <- 0
  f4 <- b1/(1 + exp(-(b4 - b2)/b3))
  # the inflection point stays at b2; the rate is b1/(4 b3) divided by 1 - f(b4)/b1
  expect_equal(unname(coef(mr)[1]), unname(b2), tolerance=1e-6)
  expect_equal(unname(coef(mr)[2]), unname(b1/(4*b3)/(1 - f4/b1)), tolerance=1e-6)
  expect_equal(unname(coef(mr)[3]), unname(1 - 0.5/(1 - f4/b1)), tolerance=1e-6)
  expect_lt(coef(mr)[[3]], 0.5)
})

test_that("bootstrap intervals of fitted groups draw their random effects at each bootstrap fit", {
  skip_on_cran()
  f <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, truncated=FALSE, data=FabaHeights25)
  mr <- maxrate(f, newdata=data.frame(Plant_ID=6))
  bt <- bootstrap(f, nsim=20, seed=1, progress=FALSE)
  ci <- suppressMessages(confint(mr, boot=bt, parm=c("tmax:Plant_ID=6", "rate:Plant_ID=6")))
  expect_true(all(is.finite(ci)) && all(ci[, 1] < ci[, 2]))
  expect_true(all(ci[, 1] <= coef(mr)[1:2] + 1 & ci[, 2] >= coef(mr)[1:2] - 1))
})
