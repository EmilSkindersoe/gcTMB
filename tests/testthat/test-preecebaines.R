# Model 1 of Preece and Baines (1978): b1 - 2 (b1 - b3) / (exp(b4 (t - b2)) + exp(b5 (t - b2))), with the
# adult height b1, the middle of puberty b2, the height b3 at b2, and the rates b4 < b5 before and
# during puberty

mpb <- 11
pbcurve <- function(t, b1, b2, b3, b4, b5) b1 - 2*(b1 - b3)/(exp(b4*(t - b2)) + exp(b5*(t - b2)))
truepb <- c(170, 13, 158, 0.1, 1)

test_that("the model index, parameters and defaults", {
  expect_equal(gcTMB:::modelindex("preecebaines"), mpb)
  expect_length(gcTMB:::unusedpars(mpb), 0)
  expect_true("preecebaines" %in% gcTMB:::modelnames)
  expect_equal(gcTMB:::modelnames[5], "beta")
  set.seed(1)
  d <- data.frame(age=rep(seq(2, 20, by=0.5), 2))
  d$height <- pbcurve(d$age, 170, 13, 158, 0.1, 1) + rnorm(nrow(d))
  f <- gcTMB(height ~ age, data=d, model="preecebaines")
  # the rates b4 and b5 are estimated by default, and the curve is not truncated
  expect_equal(f$fix, rep(FALSE, 5))
  expect_false(f$truncated)
  expect_false(any(grepl("truncated", capture.output(print(f)))))
  expect_error(gcTMB(height ~ age, data=d, model="preecebaines", b1type="rate"), "not available for the preecebaines")
})

test_that("the curve, its inverse and the peak of the pubertal growth spurt", {
  t <- c(-5, 0, 2, 8, 12, 13, 15, 20, 40, 100)
  expect_equal(gcTMB:::gcmeanR(mpb, t, 170, 13, 158, 0.1, 1), pbcurve(t, 170, 13, 158, 0.1, 1))
  q <- c(0.5, 0.8, 0.9, 0.95, 0.99, 0.999)
  tt <- gcTMB:::gctimeR(mpb, q, 13, 158, 0.1, 1, b1=170)
  expect_equal(pbcurve(tt, 170, 13, 158, 0.1, 1), q*170, tolerance=1e-10)
  # the peak height velocity: the local maximum of the growth rate after its minimum, the take-off
  for (b in list(truepb, c(180, 14, 166, 0.09, 0.8), c(165, 11, 150, 0.12, 0.9))){
    tg <- seq(b[2] - 6, b[2] + 4, by=1e-4)
    v <- diff(pbcurve(tg, b[1], b[2], b[3], b[4], b[5]))/1e-4
    # the take-off, the minimum before b2, and the peak after it
    imin <- which.min(ifelse(tg[-1] < b[2], v, Inf))
    imax <- imin - 1 + which.max(v[imin:length(v)])
    mr <- gcTMB:::maxrateR(mpb, b)
    expect_equal(unname(mr[c("tmax", "rate")]), c(tg[imax] + 5e-5, v[imax]), tolerance=1e-5)
    expect_equal(mr[["fraction"]], pbcurve(mr[["tmax"]], b[1], b[2], b[3], b[4], b[5])/b[1])
  }
  # no pubertal peak for b5/b4 below 3 + 2 sqrt(2)
  expect_true(all(is.na(gcTMB:::maxrateR(mpb, c(170, 13, 158, 0.2, 1)))))
})

test_that("Preece-Baines fits agree with nls, with timeto and maxrate", {
  set.seed(2)
  d <- data.frame(age=rep(seq(2, 20, by=0.5), 3))
  d$height <- pbcurve(d$age, 170, 13, 158, 0.1, 1) + rnorm(nrow(d))
  f <- gcTMB(height ~ age, data=d, model="preecebaines")
  expect_equal(f$opt$convergence, 0)
  n <- nls(height ~ pbcurve(age, b1, b2, b3, b4, b5), data=d, start=list(b1=170, b2=13, b3=158, b4=0.1, b5=1))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
  expect_equal(unname(predict(f, newdata=data.frame(age=c(5, 13)))), pbcurve(c(5, 13), coef(n)[1], coef(n)[2], coef(n)[3], coef(n)[4], coef(n)[5]),
               tolerance=1e-5)
  # timeto: the inverse of the curve, with delta method standard errors
  b <- bcoef(f)
  tt <- timeto(f, level=c(90, 99))
  expect_equal(unname(coef(tt)), gcTMB:::gctimeR(mpb, c(0.9, 0.99), b[2], b[3], b[4], b[5], b1=b[1]), tolerance=1e-8)
  G <- t(sapply(c(0.9, 0.99), function(q) sapply(1:5, function(j){
    h <- 1e-6*max(1, abs(b[j]))
    bp <- replace(b, j, b[j] + h); bm <- replace(b, j, b[j] - h)
    (gcTMB:::gctimeR(mpb, q, bp[2], bp[3], bp[4], bp[5], b1=bp[1]) - gcTMB:::gctimeR(mpb, q, bm[2], bm[3], bm[4], bm[5], b1=bm[1]))/(2*h)
  })))
  expect_equal(unname(sqrt(diag(vcov(tt)))), sqrt(diag(G %*% vcov(f)[1:5, 1:5] %*% t(G))), tolerance=1e-4)
  ta <- timeto(f, level=150, type="absolute")
  expect_equal(unname(pbcurve(coef(ta)[[1]], b[1], b[2], b[3], b[4], b[5])), 150, tolerance=1e-8)
  # maxrate, and profile intervals of the time and the rate of the peak
  mr <- maxrate(f)
  expect_equal(coef(mr), gcTMB:::maxrateR(mpb, b), ignore_attr=TRUE)
  ci <- suppressMessages(confint(mr, parm=c("tmax", "rate")))
  expect_true(all(ci[, 1] < coef(mr)[1:2] & coef(mr)[1:2] < ci[, 2]))
  ci <- suppressMessages(confint(tt))
  expect_true(all(ci[, 1] < coef(tt) & coef(tt) < ci[, 2]))
})
