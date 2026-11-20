# The curves that start at b4: chapmanrichards and beta (Yin et al. 2003)
#   chapmanrichards: b1 (1 - exp(-(t - b4)/b3))^c, c = exp((b2 - b4)/b3)
#   beta:            b1 (1 + (te - t)/(te - b2)) ((t - b4)/(te - b4))^((te - b4)/(te - b2)), te = b2 + b3,
#                    b1 from te on

mcr <- 8
mbeta <- 9

crcurve <- function(t, b1, b2, b3, b4) ifelse(t > b4, b1 * (1 - exp(-(t - b4)/b3))^exp((b2 - b4)/b3), 0)
betacurve <- function(t, b1, b2, b3, b4){
  te <- b2 + b3
  tc <- pmin(pmax(t, b4), te)
  ifelse(t <= b4, 0, ifelse(t >= te, b1, b1 * (1 + (te - tc)/(te - b2)) * ((tc - b4)/(te - b4))^((te - b4)/(te - b2))))
}

test_that("the model indices match the names", {
  expect_equal(gcTMB:::modelindex("chapmanrichards"), mcr)
  expect_equal(gcTMB:::modelindex("beta"), mbeta)
  expect_true(all(vapply(c(mcr, mbeta), gcTMB:::isstart, logical(1))))
  expect_length(gcTMB:::gcmodels, 9)
})

test_that("the curves are the chapmanrichards and beta curves, zero up to b4", {
  t <- c(-5, 0, 2, 5, 10, 20, 30, 45, 60, 100)
  for (b4 in c(0, 3, -4)){
    expect_equal(gcTMB:::gcmeanR(mcr, t, 50, 20, 8, b4), crcurve(t, 50, 20, 8, b4))
    expect_equal(gcTMB:::gcmeanR(mbeta, t, 50, 25, 20, b4), betacurve(t, 50, 25, 20, b4))
  }
  # truncation changes nothing
  expect_equal(gcTMB:::gcmeanR(mcr, t, 50, 20, 8, 3, trunc=TRUE), crcurve(t, 50, 20, 8, 3))
  # the beta curve stays at its asymptote from the end of growth te = b2 + b3 on
  expect_equal(gcTMB:::gcmeanR(mbeta, c(45, 60, 100), 50, 25, 20, 0), rep(50, 3))
})

test_that("chapmanrichards is the curve of growthmodels and the Richards curve with a negative asymmetry", {
  t <- seq(0.5, 80, by=0.5)
  b <- c(50, 20, 8, -2)
  cc <- exp((b[2] - b[4])/b[3])
  # growthmodels::chapmanRichards: alpha (1 - beta exp(-k t))^(1/(1 - m))
  alpha <- b[1]; k <- 1/b[3]; beta <- exp(b[4]/b[3]); m <- 1 - 1/cc
  expect_equal(gcTMB:::gcmeanR(mcr, t, b[1], b[2], b[3], b[4]), alpha * (1 - beta * exp(-k * t))^(1/(1 - m)))
  # Richards (1 + b5 exp(-z))^(-1/b5), z = (t - b2)/b3, with b5 = -1/c
  b5 <- -1/cc
  expect_equal(gcTMB:::gcmeanR(mcr, t, b[1], b[2], b[3], b[4]), b[1] * (1 + b5 * exp(-(t - b[2])/b[3]))^(-1/b5))
})

test_that("the times to reach a level are the inverse of the curves", {
  q <- c(0.01, 0.1, 0.35, 0.5, 0.9, 0.99)
  for (b4 in c(0, 3)){
    for (p in list(list(mcr, c(20, 8)), list(mcr, c(2, 8)), list(mbeta, c(25, 20)), list(mbeta, c(10, 40)))){
      tt <- gcTMB:::gctimeR(p[[1]], q, p[[2]][1] + b4, p[[2]][2], b4)
      expect_equal(gcTMB:::gcmeanR(p[[1]], tt, 1, p[[2]][1] + b4, p[[2]][2], b4), q, tolerance=1e-10)
      expect_true(all(tt > b4))
    }
  }
})

test_that("maxrate gives the inflection point, its fraction and the maximum growth rate", {
  for (p in list(list(mcr, c(50, 20, 8, 0)), list(mcr, c(50, 20, 8, 5)), list(mbeta, c(50, 25, 20, 0)), list(mbeta, c(50, 25, 10, 5)))){
    b <- p[[2]]
    t <- seq(b[4] + 1e-4, 100, by=0.001)
    f <- gcTMB:::gcmeanR(p[[1]], t, b[1], b[2], b[3], b[4])
    r <- diff(f)/0.001
    i <- which.max(r)
    mr <- gcTMB:::maxrateR(p[[1]], c(b, 1, 0))
    expect_equal(unname(mr), c(t[i] + 0.0005, r[i], mean(f[i + 0:1])/b[1]), tolerance=1e-4, label=paste(p[[1]], b[4]))
    expect_equal(mr[["tmax"]], b[2])
  }
  # closed forms: chapmanrichards (1 - 1/c)^c and b1/b3 (1 - 1/c)^(c - 1); beta with b4 = 0 the cm of
  # Yin et al. (2003)
  cc <- exp(20/8)
  expect_equal(gcTMB:::maxrateR(mcr, c(50, 20, 8, 0, 1, 0)), c(tmax=20, rate=50/8*(1 - 1/cc)^(cc - 1), fraction=(1 - 1/cc)^cc))
  te <- 45; tm <- 25
  expect_equal(gcTMB:::maxrateR(mbeta, c(50, 25, 20, 0, 1, 0))[["rate"]], 50*(2*te - tm)/(te*(te - tm))*(tm/te)^(tm/(te - tm)))
  # the additive b6 changes neither
  expect_equal(gcTMB:::maxrateR(mcr, c(50, 20, 8, 0, 1, 7)), gcTMB:::maxrateR(mcr, c(50, 20, 8, 0, 1, 0)))
})

# simulated data of a model with parameters b (b1-b4) and the additive b6
simstart <- function(mod, b, b6=0, sd=1, seed=1){
  set.seed(seed)
  d <- data.frame(day=rep(seq(1, 60, by=3), 3))
  d$y <- b6 + gcTMB:::gcmeanR(mod, d$day, b[1], b[2], b[3], b[4]) + rnorm(nrow(d), 0, sd)
  d
}

test_that("chapmanrichards and beta fits agree with nls", {
  d <- simstart(mcr, c(50, 20, 8, 0))
  f <- gcTMB(y ~ day, data=d, model="chapmanrichards")
  expect_equal(f$opt$convergence, 0)
  expect_true(f$truncated)
  expect_equal(f$fix, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  n <- nls(y ~ b1*(1 - exp(-day/b3))^exp(b2/b3), data=d, start=list(b1=50, b2=20, b3=8))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
  expect_equal(sqrt(diag(vcov(f)))[1:3], sqrt(diag(vcov(n))), tolerance=2e-2, ignore_attr=TRUE)

  d <- simstart(mbeta, c(50, 25, 20, 0))
  f <- gcTMB(y ~ day, data=d, model="beta")
  expect_equal(f$opt$convergence, 0)
  n <- nls(y ~ betacurve(day, b1, b2, b3, 0), data=d, start=list(b1=50, b2=25, b3=20))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
})

test_that("the start b4 is estimated, and the starting values are close to the estimates", {
  for (p in list(list("chapmanrichards", c(50, 23, 8, 3)), list("beta", c(50, 28, 20, 3)))){
    d <- simstart(gcTMB:::modelindex(p[[1]]), p[[2]], sd=0.5)
    f <- gcTMB(y ~ day, data=d, model=p[[1]], fix=c(b4=FALSE))
    expect_equal(f$opt$convergence, 0, label=p[[1]])
    used <- 1:4
    b <- coef(f)[paste("b", used, ":(Intercept)", sep="")]
    expect_equal(unname(b), p[[2]][used], tolerance=0.1, label=p[[1]])
    # starting values with b4 = 0 fixed, close to the estimates
    sv <- gcTMBstart(y ~ day, data=d, model=p[[1]])
    f0 <- gcTMB(y ~ day, data=d, model=p[[1]])
    u3 <- used[used < 4]
    expect_equal(unname(sv[u3]), unname(coef(f0)[paste("b", u3, ":(Intercept)", sep="")]), tolerance=0.3, label=p[[1]])
  }
})

test_that("timeto of the curves that start at b4: the inverse of the curve, with delta method standard errors", {
  for (p in list(list("chapmanrichards", c(50, 20, 8, 0)), list("beta", c(50, 25, 20, 0)))){
    mod <- gcTMB:::modelindex(p[[1]])
    d <- simstart(mod, p[[2]])
    f <- gcTMB(y ~ day, data=d, model=p[[1]])
    tt <- timeto(f, level=c(10, 50, 90))
    b <- gcTMB:::curvepars(f, f$start)
    pv <- f$start
    pv[!f$fix] <- f$estimates[seq_len(sum(!f$fix))]
    timefun <- function(par, q){
      pv[!f$fix] <- par
      gctimeR(mod, q, pv[2], pv[3], pv[4])
    }
    gctimeR <- gcTMB:::gctimeR
    est <- f$estimates[seq_len(sum(!f$fix))]
    expect_equal(unname(coef(tt)), gctimeR(mod, c(0.1, 0.5, 0.9), pv[2], pv[3], pv[4]), tolerance=1e-8, label=p[[1]])
    # delta method with numerical derivatives
    G <- t(sapply(c(0.1, 0.5, 0.9), function(q) sapply(seq_along(est), function(j){
      h <- 1e-6*max(1, abs(est[j]))
      (timefun(replace(est, j, est[j] + h), q) - timefun(replace(est, j, est[j] - h), q))/(2*h)
    })))
    V <- vcov(f)[seq_along(est), seq_along(est)]
    expect_equal(unname(sqrt(diag(vcov(tt)))), sqrt(diag(G %*% V %*% t(G))), tolerance=1e-4, label=p[[1]])
    # absolute times
    ta <- timeto(f, level=20, type="absolute")
    expect_equal(gcTMB:::gcmeanR(mod, coef(ta)[[1]], pv[1], pv[2], pv[3], pv[4]), 20, tolerance=1e-8)
    # maxrate by the delta method, and a profile interval of a time
    expect_true(all(is.finite(coef(maxrate(f))[1:2])))
  }
  d <- simstart(mcr, c(50, 20, 8, 0))
  f <- gcTMB(y ~ day, data=d, model="chapmanrichards")
  ci <- confint(timeto(f, level=50))
  expect_true(ci[1, 1] < coef(timeto(f, level=50)) && coef(timeto(f, level=50)) < ci[1, 2])
})

test_that("b1 as the maximum growth rate of the chapmanrichards and beta curves", {
  for (p in list(list("chapmanrichards", c(50, 20, 8, 0)), list("beta", c(50, 25, 20, 0)))){
    d <- simstart(gcTMB:::modelindex(p[[1]]), p[[2]])
    fa <- gcTMB(y ~ day, data=d, model=p[[1]])
    fr <- gcTMB(y ~ day, data=d, model=p[[1]], b1type="rate")
    expect_equal(logLik(fr), logLik(fa), tolerance=1e-6, label=p[[1]])
    expect_equal(coef(fr)[[1]], coef(maxrate(fa))[["rate"]], tolerance=1e-5, label=p[[1]])
  }
})
