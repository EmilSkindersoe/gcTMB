# Links: the identity link for all parameters by default, the inverse link for the scale b3, whose
# linear predictor is then the rate 1/b3, the log link for b1-b4 and the logit link for b1.

test_that("the links are the identity by default", {
  f <- gcTMB(y ~ day, truncated=FALSE, data=simgauss("logistic"))
  expect_equal(f$link, c(b1="identity", b2="identity", b3="identity", b4="identity", b5="identity", b6="identity"))
  expect_true(any(grepl("identity link for all parameters", capture.output(print(f)))))
})

test_that("the inverse link of b3 estimates the rate, with the same fit", {
  for (m in models){
    d <- simgauss(m)
    f <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m)
    fi <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m, link=c(b3="inverse"))
    expect_equal(fi$link[["b3"]], "inverse")
    expect_equal(fi$opt$objective, f$opt$objective, tolerance=1e-8, label=m)
    expect_equal(bcoef(fi)[[3]], 1/bcoef(f)[[3]], tolerance=1e-6, label=m)
    expect_equal(bcoef(fi)[1:2], bcoef(f)[1:2], tolerance=1e-6)
    # standard error of the rate by the delta method
    expect_equal(sqrt(vcov(fi)[3, 3]), sqrt(vcov(f)[3, 3])/bcoef(f)[[3]]^2, tolerance=1e-4, label=m)
    expect_equal(predict(fi), predict(f), tolerance=1e-8)
    # the parameters on their original scale and on the scale of the link
    expect_equal(unname(predict(fi, type="parameters")[1, ]), unname(bcoef(f)), tolerance=1e-6)
    expect_equal(unname(predict(fi, type="parameters", scale="link")[1, 3]), bcoef(fi)[[3]])
    s <- capture.output(print(summary(fi)))
    expect_true(any(grepl("inverse link for b3", s)))
    expect_true(any(grepl("the rate 1/b3", s)))
  }
})

test_that("the starting values of the rate are the inverse of the scale", {
  d <- simgauss("logistic")
  s <- gcTMBstart(y ~ day, truncated=FALSE, data=d)
  si <- gcTMBstart(y ~ day, truncated=FALSE, data=d, link=c(b3="inverse"))
  expect_equal(si[[3]], 1/s[[3]])
  # a given rate
  f <- gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b3="inverse"), start=c(b3=0.4), fix=c(b3=TRUE))
  expect_equal(f$start[[3]], 0.4)
  expect_equal(unname(predict(f, type="parameters")[1, "b3"]), 2.5)
})

test_that("covariates act additively on the rate with the inverse link", {
  set.seed(4)
  d <- data.frame(day=rep(days, each=4), trt=factor(rep(c("a", "b"), each=2)))
  rate <- ifelse(d$trt == "a", 0.4, 0.6)
  d$y <- 50/(1 + exp(-(d$day - 10)*rate)) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, truncated=FALSE, fform=b3 ~ trt, data=d, link=c(b3="inverse"))
  n <- nls(y ~ b1/(1 + exp(-(day - b2)*(r1 + r2*(trt == "b")))), data=d, start=list(b1=50, b2=10, r1=0.4, r2=0.2))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
})

test_that("link is checked", {
  d <- simgauss("logistic")
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link="inverse"), "named by parameter")
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b3="logg")), "Unknown link")
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b1="inverse")), "not available for b1")
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b2="logit")), "not available for b2")
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b3="logit")), "not available for b3")
  # starting values outside the range of the link
  d$day <- d$day - 30
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b2="log")), "outside the range of its log link")
  expect_error(gcTMB(y ~ day, truncated=FALSE, data=d, link=c(b7="inverse")), "parameter names")
})

test_that("log and logit links reparameterise the curve, with the same fit", {
  for (m in models){
    d <- simgauss(m)
    f <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m)
    fl <- gcTMB(y ~ day, truncated=FALSE, data=d, model=m, link=c(b1="log", b2="log", b3="log"))
    expect_equal(fl$opt$objective, f$opt$objective, tolerance=1e-8, label=m)
    expect_equal(exp(bcoef(fl)), bcoef(f), tolerance=1e-6, label=m)
    # standard errors on the log scale by the delta method
    expect_equal(sqrt(diag(vcov(fl))[1:3]), sqrt(diag(vcov(f))[1:3])/bcoef(f), tolerance=1e-4, label=m)
    expect_equal(predict(fl), predict(f), tolerance=1e-8)
    expect_equal(unname(predict(fl, type="parameters")[1, ]), unname(bcoef(f)), tolerance=1e-6)
    expect_equal(coef(timeto(fl, level=c(10, 90))), coef(timeto(f, level=c(10, 90))), tolerance=1e-6)
    s <- capture.output(print(summary(fl)))
    expect_true(any(grepl("log link for b1", s)))
    expect_true(any(grepl("b1, b2, b3 are on the log scale", s)))
  }
})

test_that("the logit link keeps a binomial asymptote in (0, 1)", {
  set.seed(4)
  d <- data.frame(day=rep(days, each=3))
  d$n <- 40
  d$s <- rbinom(nrow(d), d$n, 0.7 * plogis((d$day - 10)/2.5))
  f <- gcTMB(cbind(s, n - s) ~ day, truncated=FALSE, data=d, family="binomial")
  fl <- gcTMB(cbind(s, n - s) ~ day, truncated=FALSE, data=d, family="binomial", link=c(b1="logit"))
  expect_equal(fl$opt$objective, f$opt$objective, tolerance=1e-8)
  expect_equal(plogis(bcoef(fl)[[1]]), bcoef(f)[[1]], tolerance=1e-6)
  # confidence intervals of b1 on the logit scale, back-transformed, stay inside (0, 1)
  ci <- predict(fl, newdata=data.frame(day=0), type="parameters", interval="confidence")
  b1 <- ci[ci$parameter == "b1", ]
  expect_true(b1$lwr > 0 && b1$upr < 1)
})

test_that("covariates act multiplicatively with the log link", {
  FabaHeights25$trt <- factor(FabaHeights25$Treatment)
  f <- gcTMB(height ~ DAS, truncated=FALSE, fform=b1 ~ trt, data=FabaHeights25)
  fl <- gcTMB(height ~ DAS, truncated=FALSE, fform=b1 ~ trt, data=FabaHeights25, link=c(b1="log"))
  # a factor alone: the same fit, with ratios of the asymptotes
  expect_equal(fl$opt$objective, f$opt$objective, tolerance=1e-7)
  b <- coef(f)[1:5]
  bl <- coef(fl)[1:5]
  expect_equal(unname(exp(bl[2:5])), unname((b[1] + b[2:5])/b[1]), tolerance=1e-5)
  # random effects on the log scale of the asymptote
  fr <- gcTMB(height ~ DAS, truncated=FALSE, rform=b1 ~ 1 | plant_ID, data=FabaHeights25, link=c(b1="log"))
  expect_equal(fr$opt$convergence, 0)
  expect_true(is.finite(logLik(fr)))
})
