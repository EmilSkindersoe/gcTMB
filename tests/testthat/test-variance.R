# variance functions of the residuals: varIdent, varPower and varConstPower

d <- FabaHeights25
d$trt <- factor(d$Treatment)
d$Plant_ID <- factor(d$Plant_ID)

# exact log-likelihood of a gaussian fit without random effects and correlation
llind <- function(fit){
  mu <- predict(fit, random=FALSE)
  s <- sigma(fit) * gcTMB:::fitvarmult(fit, random=FALSE)
  sum(dnorm(fit$obj$env$data$y, mu, s, log=TRUE))
}

test_that("the constructors have the structure of those of nlme", {
  skip_if_not_installed("nlme")
  for (pair in list(list(varIdent(c(b=2), form=~ 1 | g), nlme::varIdent(c(b=2), form=~ 1 | g)),
                    list(varPower(form=~ x), nlme::varPower(form=~ x)),
                    list(varPower(0.5, form=~ x | g, fixed=c(a=1)), nlme::varPower(0.5, form=~ x | g, fixed=c(a=1))),
                    list(varConstPower(form=~ x, fixed=list(power=1)), nlme::varConstPower(form=~ x, fixed=list(power=1))))){
    expect_equal(unclass(pair[[1]]), unclass(pair[[2]]), ignore_attr="formula")
    expect_equal(class(pair[[1]]), class(pair[[2]]))
  }
  expect_error(varIdent(form=~ 1), "grouping factor")
  expect_error(varPower(form=~ 1), "covariate")
})

test_that("the likelihood is that of normal residuals with standard deviation sigma g", {
  for (v in list(varIdent(form=~ 1 | trt), varIdent(form=~ 1 | DAS), varPower(form=~ DAS), varPower(),
                 varConstPower(form=~ DAS), varConstPower(form=~ rate(.)), varPower(form=~ DAS | trt))){
    fit <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=v)
    expect_equal(fit$opt$convergence, 0)
    expect_equal(as.numeric(logLik(fit)), llind(fit), tolerance=1e-8)
    # the multiplier of the template and of R
    gv <- fit$obj$report(gcTMB:::fitpar(fit))$gv
    expect_equal(as.vector(gv), gcTMB:::fitvarmult(fit), tolerance=1e-6)
  }
})

test_that("the fits agree with nlme::gnls", {
  skip_if_not_installed("nlme")
  for (k in 1:4){
    vg <- list(varIdent(form=~ 1 | trt), varPower(form=~ DAS), varConstPower(form=~ DAS), varIdent(form=~ 1 | DAS))[[k]]
    vn <- list(nlme::varIdent(form=~ 1 | trt), nlme::varPower(form=~ DAS), nlme::varConstPower(form=~ DAS),
               nlme::varIdent(form=~ 1 | DAS))[[k]]
    fit <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=vg)
    cf <- coef(fit)
    g <- nlme::gnls(Height ~ b1/(1 + exp(-(DAS - b2)/b3)), data=d, start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]), weights=vn)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(g)), tolerance=1e-6)
  }
  # nlme's objects are accepted
  f1 <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=nlme::varPower(form=~ DAS))
  f2 <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=varPower(form=~ DAS))
  expect_equal(logLik(f1), logLik(f2))
  # with a correlation structure
  fs <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=varPower(form=~ DAS), correlation=corCAR1(form=~ DAS | Plant_ID))
  cf <- coef(fs)
  gs <- nlme::gnls(Height ~ b1/(1 + exp(-(DAS - b2)/b3)), data=d, start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]),
                   weights=nlme::varPower(form=~ DAS), correlation=nlme::corCAR1(form=~ DAS | Plant_ID))
  expect_equal(as.numeric(logLik(fs)), as.numeric(logLik(gs)), tolerance=1e-5)
})

test_that("varIdent ratios are relative to the first level and fixed values stay fixed", {
  fit <- gcTMB(Height ~ DAS, data=d, variance=varIdent(form=~ 1 | trt))
  vp <- varpars(fit)
  expect_equal(names(vp), levels(d$trt))
  expect_equal(vp[[1]], 1)
  fx <- gcTMB(Height ~ DAS, data=d, variance=varPower(form=~ DAS, fixed=1))
  expect_equal(varpars(fx), c(power=1))
  expect_false(any(grepl("var_par", names(coef(fx)))))
  fc <- gcTMB(Height ~ DAS, data=d, variance=varConstPower(form=~ DAS, fixed=list(power=1)))
  expect_equal(varpars(fc)[["power"]], 1)
  expect_true("var_par:const" %in% names(coef(fc)))
})

test_that("residuals, simulation, prediction and update use the variance function", {
  fit <- gcTMB(Height ~ DAS, data=d, variance=varPower(form=~ DAS))
  g <- d$DAS^varpars(fit)[["power"]]
  r <- residuals(fit, type="pearson")
  expect_equal(r, (d$Height - fitted(fit))/(sigma(fit)*g))
  pr <- predict(fit, newdata=data.frame(DAS=c(20, 60)), interval="prediction")
  se <- predict(fit, newdata=data.frame(DAS=c(20, 60)), se.fit=TRUE)$se.fit
  sdp <- sqrt(se^2 + (sigma(fit)*c(20, 60)^varpars(fit)[["power"]])^2)
  expect_equal(pr$upr - pr$fit, qnorm(0.975)*sdp, tolerance=1e-8)
  # simulated residuals have the standard deviation of the variance function
  s <- simulate(fit, nsim=2000, seed=1, newdata=data.frame(DAS=c(14, 61)))
  sdsim <- apply(as.matrix(s) - as.matrix(attr(s, "mean")), 1, sd)
  expect_equal(sdsim, sigma(fit)*c(14, 61)^varpars(fit)[["power"]], tolerance=0.05)
  s2 <- simulate(fit, nsim=2000, seed=1)
  sd2 <- apply(as.matrix(s2) - as.matrix(attr(s2, "mean")), 1, sd)
  expect_equal(sd2, sigma(fit)*g, tolerance=0.1)
  expect_equal(varpars(update(fit, varpar=c(power=1))), c(power=1))
  expect_error(update(fit, varpar=c(const=1)), "varpar must be")
  ft <- gcTMB(Height ~ DAS, data=d, variance=varIdent(form=~ 1 | DAS))
  expect_error(simulate(ft, newdata=data.frame(DAS=100)), "no variance for the strata")
  expect_error(gcTMB(Pods ~ DAS, data=FabaCounts23, family="poisson", variance=varPower()), "gaussian and lognormal")
})

test_that("timeto with residual=TRUE uses the variance function at the estimated time", {
  for (v in list(varPower(form=~ DAS), varPower(), varConstPower(form=~ rate(.)))){
    fit <- gcTMB(Height ~ DAS, data=d, variance=v)
    t0 <- timeto(fit, level=c(10, 50, 90))
    t1 <- timeto(fit, level=c(10, 50, 90), residual=TRUE)
    tt <- coef(t0)
    pd <- gcTMB:::preddata(fit, data.frame(DAS=tt), random=FALSE)
    h <- 1e-5*(1 + abs(tt))
    p1 <- p2 <- pd
    p1$x <- tt + h
    p2$x <- tt - h
    fr <- (gcTMB:::predmean(p1, pd$eta) - gcTMB:::predmean(p2, pd$eta))/(2*h)
    g <- gcTMB:::varmult(fit$variance, fit$obj$env$parList(fit$opt$par)$var_par, fit$varlevels, data.frame(DAS=tt),
                         gcTMB:::predmean(pd, pd$eta), function() fr)
    expect_equal(unname(diag(vcov(t1)) - diag(vcov(t0))), unname((sigma(fit)*g/fr)^2), tolerance=1e-6)
  }
  # a single curve: the times and growth rates are the same, the residual variances differ by g^2
  fi <- gcTMB(Height ~ DAS, data=d, variance=varIdent(form=~ 1 | trt))
  nd <- data.frame(trt=levels(d$trt))
  dv <- diag(vcov(timeto(fi, newdata=nd, residual=TRUE))) - diag(vcov(timeto(fi, newdata=nd)))
  expect_equal(unname(dv/dv[1]), unname(varpars(fi)^2), tolerance=1e-8)
  ft <- gcTMB(Height ~ DAS, data=d, variance=varIdent(form=~ 1 | DAS))
  expect_error(timeto(ft, residual=TRUE), "cannot be evaluated")
})

test_that("varExp and varComb have the structure of those of nlme", {
  skip_if_not_installed("nlme")
  for (pair in list(list(varExp(form=~ x), nlme::varExp(form=~ x)),
                    list(varExp(0.1, form=~ x | g, fixed=c(a=0)), nlme::varExp(0.1, form=~ x | g, fixed=c(a=0))))){
    expect_equal(unclass(pair[[1]]), unclass(pair[[2]]), ignore_attr="formula")
    expect_equal(class(pair[[1]]), class(pair[[2]]))
  }
  vc <- varComb(varIdent(form=~ 1 | trt), varPower())
  expect_equal(class(vc), class(nlme::varComb(nlme::varIdent(form=~ 1 | trt), nlme::varPower())))
  expect_equal(names(vc), c("A", "B"))
  # nested combinations are flattened
  expect_equal(length(varComb(vc, varExp())), 3)
  expect_error(varComb(varPower()), "at least two")
  expect_error(varPoly(degree=0), "positive whole number")
  expect_error(varPoly(c(1, 2, 3)), "one value for each power")
})

test_that("varExp, varPoly and varComb have the likelihood of sigma g, with g the product of the components", {
  vs <- list(varExp(), varExp(form=~ DAS), varExp(form=~ DAS | trt), varPoly(), varPoly(degree=1, form=~ DAS),
             varPoly(form=~ DAS | trt), varPoly(form=~ rate(.)), varComb(varPower(), varExp()),
             varComb(varIdent(form=~ 1 | trt), varPoly(form=~ DAS)))
  for (v in vs){
    fit <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=v)
    expect_equal(fit$opt$convergence, 0, label=deparse1(v))
    expect_equal(as.numeric(logLik(fit)), llind(fit), tolerance=1e-8)
    gv <- fit$obj$report(gcTMB:::fitpar(fit))$gv
    expect_equal(as.vector(gv), gcTMB:::fitvarmult(fit), tolerance=1e-6)
  }
  # the product: g of varComb(varPower(), varExp()) is |f|^delta exp(t f)
  fc <- gcTMB(Height ~ DAS, data=d, variance=varComb(varPower(), varExp()))
  vp <- varpars(fc)
  expect_equal(names(vp), c("A.power", "B.expon"))
  f <- fitted(fc)
  expect_equal(gcTMB:::fitvarmult(fc), abs(f)^vp[["A.power"]]*exp(vp[["B.expon"]]*f), tolerance=1e-8)
  expect_match(paste(capture.output(print(fc)), collapse=" "), "product of A. power of fitted")
})

test_that("varExp and varComb agree with nlme::gnls", {
  skip_if_not_installed("nlme")
  vg <- list(varExp(form=~ DAS), varComb(varIdent(form=~ 1 | trt), varExp(form=~ DAS)), varComb(varPower(form=~ DAS), varExp(form=~ DAS)))
  vn <- list(nlme::varExp(form=~ DAS), nlme::varComb(nlme::varIdent(form=~ 1 | trt), nlme::varExp(form=~ DAS)),
             nlme::varComb(nlme::varPower(form=~ DAS), nlme::varExp(form=~ DAS)))
  for (k in seq_along(vg)){
    fit <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=vg[[k]])
    cf <- coef(fit)
    g <- nlme::gnls(Height ~ b1/(1 + exp(-(DAS - b2)/b3)), data=d, start=c(b1=cf[[1]], b2=cf[[2]], b3=cf[[3]]), weights=vn[[k]])
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(g)), tolerance=1e-6)
  }
  # nlme's varComb is accepted
  f1 <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=vn[[2]])
  f2 <- gcTMB(Height ~ DAS, data=d, truncated=FALSE, variance=vg[[2]])
  expect_equal(logLik(f1), logLik(f2))
})

test_that("varPoly: a polynomial standard deviation centred at the mean, with a barrier where it is 0 or below", {
  fit <- gcTMB(Height ~ DAS, data=d, variance=varPoly(form=~ DAS))
  vp <- varpars(fit)
  expect_equal(names(vp), c("poly1", "poly2"))
  m <- mean(d$DAS)
  expect_equal(gcTMB:::fitvarmult(fit), 1 + vp[["poly1"]]*(d$DAS - m) + vp[["poly2"]]*(d$DAS - m)^2, tolerance=1e-8)
  expect_match(paste(capture.output(print(fit)), collapse=" "), "polynomial of degree 2 in DAS - ")
  # the same fit from other starting values; fixed coefficients stay fixed
  f2 <- gcTMB(Height ~ DAS, data=d, variance=varPoly(c(0.01, 0), form=~ DAS))
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(fit)), tolerance=1e-6)
  f1 <- gcTMB(Height ~ DAS, data=d, variance=varPoly(form=~ DAS, fixed=c(poly2=0)))
  expect_equal(varpars(f1)[["poly2"]], 0)
  expect_equal(attr(logLik(f1), "df"), attr(logLik(fit), "df") - 1)
  # Wald intervals on the scale of the coefficients
  ci <- confint(fit)
  expect_true(all(c("poly1", "poly2") %in% rownames(ci)))
  expect_true(ci["poly1", 1] < vp[["poly1"]] && vp[["poly1"]] < ci["poly1", 2])
  # the barrier: a polynomial 0 or below at an observation gives a log-likelihood of minus infinity
  fb <- update(fit, varpar=c(poly2=-1))
  expect_equal(as.numeric(logLik(fb)), -Inf)
  expect_error(predict(fb, interval="prediction"), "0 or below")
  # update and refit with the estimates
  expect_equal(varpars(update(fit, varpar=c(poly1=0.02)))[["poly1"]], 0.02)
  fr <- refit(fit, startfrom="fit")
  expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fit)), tolerance=1e-8)
  expect_equal(eval(fr$call$variance)[["poly1"]], vp[["poly1"]], ignore_attr=TRUE)
  # by a grouping factor
  fg <- gcTMB(Height ~ DAS, data=d, variance=varPoly(degree=1, form=~ DAS | trt))
  expect_equal(names(varpars(fg)), paste("poly1:", levels(d$trt), sep=""))
})

test_that("predict without newdata works with integer times", {
  fit <- gcTMB(Height ~ DAS, data=d, variance=varPower(form=~ DAS))
  expect_true(is.integer(d$DAS))
  p <- predict(fit, interval="prediction")
  expect_equal(nrow(p), nrow(d))
  expect_equal(p$fit, predict(fit, newdata=transform(d, DAS=as.numeric(DAS)), interval="prediction")$fit)
})
