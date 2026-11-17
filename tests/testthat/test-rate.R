# b1 as the maximum growth rate (b1type = "rate"): the same models as with the asymptote, with b1 the
# rate of maxrate()

test_that("b1type = 'rate' gives the same fit with b1 the maximum growth rate", {
  cases <- list(list(model="logistic"), list(model="gompertz"), list(model="hill"), list(model="weibull"),
                list(model="normal"), list(model="gamma"), list(model="logistic", truncated=FALSE),
                list(model="logistic", fix=c(b5=FALSE), link=c(b5="log")),
                list(model="hill", fix=c(b5=FALSE), link=c(b5="log")),
                list(model="logistic", fix=c(b4=FALSE)))
  for (cs in cases){
    lab <- paste(deparse(cs), collapse="")
    fa <- do.call(gcTMB, c(list(height ~ DAS, data=FabaHeights25), cs))
    fr <- do.call(gcTMB, c(list(height ~ DAS, data=FabaHeights25, b1type="rate"), cs))
    expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fa)), tolerance=1e-6, label=lab)
    expect_equal(bcoef(fr)[[1]], coef(maxrate(fa))[["rate"]], tolerance=1e-4, label=lab)
    # the same curve, and maxrate of the rate fit is its b1
    expect_equal(fitted(fr), fitted(fa), tolerance=1e-4, label=lab)
    expect_equal(coef(maxrate(fr))[["rate"]], bcoef(fr)[[1]], tolerance=1e-6, label=lab)
  }
})

test_that("b1 as the rate: predictions, times, parameters and summary", {
  d <- FabaHeights25
  d$trt <- factor(d$Treatment)
  fa <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=d)
  fr <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=d, b1type="rate")
  nd <- data.frame(trt=levels(d$trt), DAS=30)
  expect_equal(predict(fr, newdata=nd), predict(fa, newdata=nd), tolerance=1e-5)
  expect_equal(predict(fr, newdata=nd, se.fit=TRUE)$se.fit, predict(fa, newdata=nd, se.fit=TRUE)$se.fit, tolerance=1e-3)
  for (type in c("relative", "absolute")){
    ta <- timeto(fa, level=50, newdata=nd["trt"], type=type)
    tr <- timeto(fr, level=50, newdata=nd["trt"], type=type)
    expect_equal(coef(tr), coef(ta), tolerance=1e-4, label=type)
    expect_equal(sqrt(diag(vcov(tr))), sqrt(diag(vcov(ta))), tolerance=1e-2, label=type)
  }
  pa <- predict(fa, newdata=nd, type="parameters")
  pr <- predict(fr, newdata=nd, type="parameters")
  expect_equal(unname(pr[, "asymptote"]), unname(pa[, "b1"]), tolerance=1e-4)
  expect_equal(unname(pr[, "b1"]), unname(coef(maxrate(fa, newdata=nd["trt"]))[paste("rate:trt=", levels(d$trt), sep="")]), tolerance=1e-4)
  expect_output(print(summary(fr)), "b1 is the maximum growth rate")
  # profile of the rate: the profile of the coefficient b1
  ci <- suppressMessages(confint(maxrate(fr, newdata=nd["trt"]), parm="rate:trt=0.1"))
  expect_equal(unname(ci[1, ]), unname(confint(fr, parm=1)[1, ]), tolerance=1e-3)
  expect_error(confint(timeto(fr, level=50, newdata=nd["trt"][1, , drop=FALSE], type="absolute")), "not available with b1type")
})

test_that("b1 as the rate needs b5 = 1 for the gompertz, weibull, normal and gamma models", {
  for (m in c("gompertz", "weibull", "normal", "gamma")){
    expect_error(gcTMB(height ~ DAS, data=FabaHeights25, model=m, b1type="rate", fix=c(b5=FALSE)), "b5 must be fixed at 1", label=m)
    expect_error(gcTMB(height ~ DAS, data=FabaHeights25, model=m, b1type="rate", start=c(b5=2)), "b5 must be fixed at 1", label=m)
  }
  expect_silent(gcTMB(height ~ DAS, data=FabaHeights25, model="gompertz", b1type="rate", start=c(b5=1)))
})

test_that("pairwise compares a parameter between curves", {
  d <- FabaHeights25
  d$trt <- factor(d$Treatment)
  fr <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=d, b1type="rate", link=c(b1="log"))
  nd <- data.frame(trt=levels(d$trt))
  pw <- pairwise(fr, nd, adjust="none")
  expect_s3_class(pw, "pairwise.gcTMB")
  expect_equal(nrow(pw), choose(5, 2))
  expect_equal(rownames(pw)[1], "trt=0.1 - trt=0.3")
  # against the reference level: the coefficients of the treatments, with their standard errors
  cf <- coef(fr)
  se <- sqrt(diag(vcov(fr)))
  expect_equal(pw$Estimate[1:4], -unname(cf[paste("b1:trt", levels(d$trt)[-1], sep="")]))
  expect_equal(pw$Std.Error[1:4], unname(se[paste("b1:trt", levels(d$trt)[-1], sep="")]), tolerance=1e-6)
  expect_equal(pw$p.value[1], 2*pt(-abs(pw$Estimate[1]/pw$Std.Error[1]), fr$rdf))
  expect_equal(pw$ratio, exp(pw$Estimate))
  # adjusted p-values
  expect_equal(pairwise(fr, nd, adjust="holm")$p.value, p.adjust(pw$p.value, "holm"))
  expect_output(print(pw), "maximum growth rate")
  expect_output(print(pw), "trt=0.1 - trt=0.3")
  expect_error(pairwise(fr, nd[1, , drop=FALSE]), "at least two rows")
  expect_error(pairwise(fr, nd, parm="b3"), "the same for all rows")
})

test_that("the single-step (Tukey) adjustment", {
  # balanced groups with independent means of equal variance: the largest absolute pairwise
  # statistic is the studentized range divided by sqrt(2)
  k <- 4
  pr <- utils::combn(k, 2)
  D <- t(apply(pr, 2, function(ij) replace(numeric(k), ij, c(1, -1))))
  R <- cov2cor(D %*% t(D))
  for (df in list(NULL, 15)){
    mx <- gcTMB:::maxabsdraws(R, df, nsim=2e5)
    for (cc in c(1.5, 2.5, 3)){
      expect_equal(mean(mx <= cc), ptukey(cc*sqrt(2), k, if (is.null(df)) Inf else df), tolerance=0.005)
    }
  }

  d <- FabaHeights25
  d$trt <- factor(d$Treatment)
  fr <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=d, b1type="rate", link=c(b1="log"))
  nd <- data.frame(trt=levels(d$trt))
  set.seed(5)
  r1 <- runif(1)
  set.seed(5)
  pt <- pairwise(fr, nd)
  expect_equal(runif(1), r1)
  expect_equal(pairwise(fr, nd)$p.value, pt$p.value)
  p0 <- pairwise(fr, nd, adjust="none")
  pb <- pairwise(fr, nd, adjust="bonferroni")
  expect_true(all(pt$p.value >= p0$p.value - 0.002 & pt$p.value <= pb$p.value + 0.002))
  # simultaneous intervals: wider than the unadjusted ones, narrower than Bonferroni
  crit <- attr(pt, "critical")
  expect_true(crit > qt(0.975, fr$rdf) && crit < qt(1 - 0.025/nrow(pt), fr$rdf))
  expect_equal(pt$upper, pt$Estimate + crit*pt$Std.Error)
  expect_output(print(pt), "single-step")
})
