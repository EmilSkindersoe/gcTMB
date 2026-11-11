# Links: the identity link for all parameters by default, and the inverse link for the scale b3,
# whose linear predictor is then the rate 1/b3.

test_that("the links are the identity by default", {
  f <- gcTMB(y ~ day, data=simgauss("logistic"))
  expect_equal(f$link, c(b1="identity", b2="identity", b3="identity", b4="identity", b5="identity"))
  expect_true(any(grepl("identity link for all parameters", capture.output(print(f)))))
})

test_that("the inverse link of b3 estimates the rate, with the same fit", {
  for (m in models){
    d <- simgauss(m)
    f <- gcTMB(y ~ day, data=d, model=m)
    fi <- gcTMB(y ~ day, data=d, model=m, link=c(b3="inverse"))
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
  s <- gcTMBstart(y ~ day, data=d)
  si <- gcTMBstart(y ~ day, data=d, link=c(b3="inverse"))
  expect_equal(si[[3]], 1/s[[3]])
  # a given rate
  f <- gcTMB(y ~ day, data=d, link=c(b3="inverse"), start=c(b3=0.4), fix=c(b3=TRUE))
  expect_equal(f$start[[3]], 0.4)
  expect_equal(unname(predict(f, type="parameters")[1, "b3"]), 2.5)
})

test_that("covariates act additively on the rate with the inverse link", {
  set.seed(4)
  d <- data.frame(day=rep(days, each=4), trt=factor(rep(c("a", "b"), each=2)))
  rate <- ifelse(d$trt == "a", 0.4, 0.6)
  d$y <- 50/(1 + exp(-(d$day - 10)*rate)) + rnorm(nrow(d), 0, 1)
  f <- gcTMB(y ~ day, fform=b3 ~ trt, data=d, link=c(b3="inverse"))
  n <- nls(y ~ b1/(1 + exp(-(day - b2)*(r1 + r2*(trt == "b")))), data=d, start=list(b1=50, b2=10, r1=0.4, r2=0.2))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
})

test_that("link is checked", {
  d <- simgauss("logistic")
  expect_error(gcTMB(y ~ day, data=d, link="inverse"), "named by parameter")
  expect_error(gcTMB(y ~ day, data=d, link=c(b3="logg")), "Unknown link")
  expect_error(gcTMB(y ~ day, data=d, link=c(b1="inverse")), "not available for b1")
  expect_error(gcTMB(y ~ day, data=d, link=c(b2="log")), "Unknown link")
  expect_error(gcTMB(y ~ day, data=d, link=c(b6="inverse")), "parameter names")
})
