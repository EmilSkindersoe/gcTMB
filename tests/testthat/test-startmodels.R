# The chapmanrichards, beta (Yin et al. 2003) and korf (Kuzelka and Marusak 2015) curves, never truncated:
#   chapmanrichards: b4 + b1 (1 - exp(-b2 t))^b3 (Richards 1959), with the rate b2, the shape b3 and
#                    the intercept b4
#   beta_capped:     b1 (1 + (b3 - t)/(b3 - b2)) ((t - b4)/(b3 - b4))^((b3 - b4)/(b3 - b2)), with the time
#                    b2 = tm of the maximum growth rate, the end of growth b3 = te and the start b4 = tb;
#                    zero up to b4 and b1 from b3 on
#   korf:            b4 + b1 exp(b2/((1 - b3) t^(b3 - 1))), as in Kuzelka and Marusak (2015), with the
#                    asymptote b1, the coefficient b2 > 0, the shape b3 > 1 and the intercept b4

mcr <- 8
mbeta <- 9
mko <- 10

crcurve <- function(t, b1, b2, b3, b4=0) b4 + ifelse(t > 0, b1 * (1 - exp(-b2*t))^b3, 0)
kocurve <- function(t, b1, b2, b3, b4=0) b4 + ifelse(t > 0, b1 * exp(b2/((1 - b3)*pmax(t, 1e-300)^(b3 - 1))), 0)
betacurve <- function(t, b1, b2, b3, b4){
  te <- b3
  tc <- pmin(pmax(t, b4), te)
  ifelse(t <= b4, 0, ifelse(t >= te, b1, b1 * (1 + (te - tc)/(te - b2)) * ((tc - b4)/(te - b4))^((te - b4)/(te - b2))))
}

test_that("the model indices match the names", {
  expect_equal(gcTMB:::modelindex("chapmanrichards"), mcr)
  expect_equal(gcTMB:::modelindex("beta_capped"), mbeta)
  expect_equal(gcTMB:::modelindex("beta"), 12)
  expect_equal(gcTMB:::modelindex("korf"), mko)
  expect_true(all(vapply(c(mcr, mbeta, mko), gcTMB:::isstart, logical(1))))
  expect_equal(gcTMB:::unusedpars(mcr), 5)
  expect_equal(gcTMB:::unusedpars(mko), 5)
  expect_equal(gcTMB:::unusedpars(mbeta), 4:5)
})

test_that("the curves are the chapmanrichards, beta and korf curves", {
  t <- c(-5, -0.5, 0, 2, 5, 10, 20, 30, 45, 60, 100)
  for (b4 in c(0, 1.3)){
    expect_equal(gcTMB:::gcmeanR(mcr, t, 50, 0.06, 3, b4), crcurve(t, 50, 0.06, 3, b4))
    expect_equal(gcTMB:::gcmeanR(mko, t, 50, 15, 2, b4), kocurve(t, 50, 15, 2, b4))
    expect_equal(gcTMB:::gcmeanR(mko, t, 50, 3.6, 1.6, b4), kocurve(t, 50, 3.6, 1.6, b4))
  }
  for (b4 in c(0, 3, -4)){
    expect_equal(gcTMB:::gcmeanR(mbeta, t, 50, 25, 45, b4), betacurve(t, 50, 25, 45, b4))
  }
  # truncation changes nothing
  expect_equal(gcTMB:::gcmeanR(mcr, t, 50, 0.06, 3, 1.3, trunc=TRUE), crcurve(t, 50, 0.06, 3, 1.3))
  # the beta curve stays at its asymptote from the end of growth b3 on
  expect_equal(gcTMB:::gcmeanR(mbeta, c(45, 60, 100), 50, 25, 45, 0), rep(50, 3))
  # with b4 = 0, eqn (8) of Yin et al. (2003)
  tt <- seq(1, 44, by=1)
  expect_equal(gcTMB:::gcmeanR(mbeta, tt, 50, 25, 45, 0), 50*(1 + (45 - tt)/(45 - 25))*(tt/45)^(45/(45 - 25)))
})

test_that("chapmanrichards is model M8 of Ciceu et al. (2020), and korf the Korf curve", {
  D <- seq(1, 80, by=0.5)
  # M8: 1.3 + a (1 - exp(-k D))^c
  a <- 30; k <- 0.05; cc <- 3
  expect_equal(gcTMB:::gcmeanR(mcr, D, a, k, cc, 1.3), 1.3 + a*(1 - exp(-k*D))^cc)
  # Korf, as in Kuzelka and Marusak (2015): 1.3 + a exp(b/((1 - c) D^(c - 1)))
  a <- 30; b <- 8; cc <- 1.7
  expect_equal(gcTMB:::gcmeanR(mko, D, a, b, cc, 1.3), 1.3 + a*exp(b/((1 - cc)*D^(cc - 1))))
})

test_that("the times to reach a level are the inverse of the curves", {
  q <- c(0.01, 0.1, 0.35, 0.5, 0.9, 0.99)
  for (b4 in c(0, 1.3)){
    for (p in list(list(mcr, c(0.06, 3)), list(mcr, c(0.2, 0.8)), list(mko, c(15, 2)), list(mko, c(0.5, 2)), list(mko, c(3.6, 1.6)))){
      tt <- gcTMB:::gctimeR(p[[1]], q, p[[2]][1], p[[2]][2], b4)
      expect_equal(gcTMB:::gcmeanR(p[[1]], tt, 1, p[[2]][1], p[[2]][2], b4) - b4, q, tolerance=1e-10)
    }
  }
  for (b4 in c(0, 3)){
    for (b in list(c(25, 20), c(10, 40))){
      tt <- gcTMB:::gctimeR(mbeta, q, b[1] + b4, b[1] + b[2] + b4, b4)
      expect_equal(gcTMB:::gcmeanR(mbeta, tt, 1, b[1] + b4, b[1] + b[2] + b4, b4), q, tolerance=1e-10)
      expect_true(all(tt > b4))
    }
  }
  # b2 from a time and from the inflection point, for the profiles
  for (p in list(list(mcr, c(0.06, 3)), list(mko, c(15, 2)), list(mko, c(3.6, 1.6)))){
    m <- gcTMB:::gcmodel(p[[1]])
    t50 <- gcTMB:::gctimeR(p[[1]], 0.5, p[[2]][1], p[[2]][2], 0)
    expect_equal(m$b2fromt(0.5, t50, p[[2]][2]), p[[2]][1])
    tm <- gcTMB:::maxrateR(p[[1]], c(1, p[[2]], 0, 1))[["tmax"]]
    expect_equal(m$b2fromtmax(tm, p[[2]][2]), p[[2]][1])
  }
})

test_that("maxrate gives the inflection point, its fraction and the maximum growth rate", {
  for (p in list(list(mcr, c(50, 0.06, 3, 0)), list(mcr, c(50, 0.06, 3, 1.3)), list(mko, c(50, 22, 2, 0)), list(mko, c(50, 8, 2, 1.3)),
                 list(mko, c(50, 3.6, 1.6, 0)), list(mbeta, c(50, 25, 45, 0)), list(mbeta, c(50, 25, 35, 5)))){
    b <- p[[2]]
    t <- seq(max(b[4], 0) + 1e-4, 100, by=0.001)
    f <- gcTMB:::gcmeanR(p[[1]], t, b[1], b[2], b[3], b[4]) - gcTMB:::intercept(p[[1]], b[4])
    r <- diff(f)/0.001
    i <- which.max(r)
    mr <- gcTMB:::maxrateR(p[[1]], c(b, 1))
    expect_equal(unname(mr), c(t[i] + 0.0005, r[i], mean(f[i + 0:1])/b[1]), tolerance=1e-4, label=paste(p[[1]], b[4]))
  }
  # closed forms: chapmanrichards at log(c)/k, (1 - 1/c)^c and b1 k (1 - 1/c)^(c - 1); korf with b3 = 2
  # at b2/2, exp(-2) and 4 b1 exp(-2)/b2; beta with b4 = 0 the cm of Yin et al. (2003)
  expect_equal(gcTMB:::maxrateR(mcr, c(50, 0.06, 3, 0, 1)), c(tmax=log(3)/0.06, rate=50*0.06*(2/3)^2, fraction=(2/3)^3))
  expect_equal(gcTMB:::maxrateR(mko, c(50, 22, 2, 0, 1)), c(tmax=11, rate=4*50*exp(-2)/22, fraction=exp(-2)))
  te <- 45; tm <- 25
  expect_equal(gcTMB:::maxrateR(mbeta, c(50, 25, 45, 0, 1))[["rate"]], 50*(2*te - tm)/(te*(te - tm))*(tm/te)^(tm/(te - tm)))
  # the intercept changes neither
  expect_equal(gcTMB:::maxrateR(mcr, c(50, 0.06, 3, 1.3, 1)), gcTMB:::maxrateR(mcr, c(50, 0.06, 3, 0, 1)))
  # chapmanrichards without an inflection point at a positive diameter (b3 <= 1)
  expect_equal(gcTMB:::maxrateR(mcr, c(50, 0.06, 1, 0, 1)), c(tmax=0, rate=50*0.06, fraction=0))
})

# simulated data of a model with parameters b (b1-b4, b1 the asymptote)
simstart <- function(mod, b, sd=1, seed=1){
  set.seed(seed)
  d <- data.frame(day=rep(seq(1, 60, by=3), 3))
  d$y <- gcTMB:::gcmeanR(mod, d$day, b[1], b[2], b[3], b[4]) + rnorm(nrow(d), 0, sd)
  d
}

test_that("chapmanrichards, beta and korf fits agree with nls", {
  d <- simstart(mcr, c(50, 0.06, 3, 0))
  f <- gcTMB(y ~ day, data=d, model="chapmanrichards", start=c(b4=0))
  expect_equal(f$opt$convergence, 0)
  expect_false(f$truncated)
  expect_equal(f$fix, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  n <- nls(y ~ b1*(1 - exp(-b2*day))^b3, data=d, start=list(b1=50, b2=0.06, b3=3))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
  expect_equal(sqrt(diag(vcov(f)))[1:3], sqrt(diag(vcov(n))), tolerance=2e-2, ignore_attr=TRUE)

  # korf: the asymptote b1, the coefficient b2 and the shape b3 estimated
  d <- simstart(mko, c(50, 3.6, 1.6, 0), sd=0.5)
  f <- gcTMB(y ~ day, data=d, model="korf", start=c(b4=0))
  expect_equal(f$opt$convergence, 0)
  expect_false(f$truncated)
  expect_equal(f$fix, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(f$link[["b1"]], "identity")
  n <- nls(y ~ b1*exp(b2/((1 - b3)*day^(b3 - 1))), data=d, start=list(b1=50, b2=3.6, b3=1.6))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-4)
  # the shape fixed at 2
  d1 <- simstart(mko, c(50, 22, 2, 0))
  f1 <- gcTMB(y ~ day, data=d1, model="korf", start=c(b3=2, b4=0), fix=c(b3=TRUE))
  expect_equal(f1$fix, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  n1 <- nls(y ~ b1*exp(-b2/day), data=d1, start=list(b1=50, b2=22))
  expect_equal(unname(bcoef(f1)), unname(coef(n1)), tolerance=1e-5)

  d <- simstart(mbeta, c(50, 25, 45, 0))
  f <- gcTMB(y ~ day, data=d, model="beta_capped")
  expect_equal(f$opt$convergence, 0)
  n <- nls(y ~ betacurve(day, b1, b2, b3, 0), data=d, start=list(b1=50, b2=25, b3=45))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
})

test_that("the intercept b4 of chapmanrichards and korf, fixed and estimated", {
  set.seed(4)
  d <- data.frame(dbh=runif(150, 5, 60))
  d$h <- gcTMB:::gcmeanR(mcr, d$dbh, 30, 0.06, 2, 1.3) + rnorm(150, 0, 1)
  # fixed at the breast height
  f <- gcTMB(h ~ dbh, data=d, model="chapmanrichards", start=c(b4=1.3))
  n <- nls(h ~ 1.3 + b1*(1 - exp(-b2*dbh))^b3, data=d, start=list(b1=30, b2=0.06, b3=2))
  expect_equal(unname(bcoef(f)[1:3]), unname(coef(n)), tolerance=1e-5)
  expect_equal(unname(predict(f, newdata=data.frame(dbh=0))), 1.3)
  expect_false(any(grepl("truncated", capture.output(print(f)))))
  # estimated
  fe <- gcTMB(h ~ dbh, data=d, model="chapmanrichards", fix=c(b4=FALSE))
  ne <- nls(h ~ b4 + b1*(1 - exp(-b2*dbh))^b3, data=d, start=list(b1=30, b2=0.06, b3=2, b4=1.3))
  expect_equal(unname(bcoef(fe)), unname(coef(ne)), tolerance=1e-4)
  # korf from the breast height
  d$h <- gcTMB:::gcmeanR(mko, d$dbh, 30, 14, 2, 1.3) + rnorm(150, 0, 1)
  fw <- gcTMB(h ~ dbh, data=d, model="korf", start=c(b4=1.3))
  nw <- nls(h ~ 1.3 + a*exp(b/((1 - c)*dbh^(c - 1))), data=d, start=list(a=30, b=14, c=2))
  expect_equal(unname(bcoef(fw)), unname(coef(nw)), tolerance=1e-4)
  # timeto and maxrate do not depend on the intercept; absolute times include it
  cw <- coef(fw)
  expect_equal(coef(timeto(fw, level=50)), unname(gcTMB:::gctimeR(mko, 0.5, cw[[2]], cw[[3]], 1.3)), ignore_attr=TRUE)
  ta <- timeto(fw, level=10, type="absolute")
  expect_equal(gcTMB:::gcmeanR(mko, coef(ta)[[1]], cw[[1]], cw[[2]], cw[[3]], 1.3), 10, tolerance=1e-8)
  p <- plot(fw, which=2, draw=FALSE)
  expect_equal(p$time$labels$title, "Residuals vs dbh")
})

test_that("b4 defaults to the breast height for chapmanrichards and korf, to 0 for the other curves, and is estimated for preecebaines", {
  f <- gcTMB(height ~ dbh, data=PineHeightDiameter, model="chapmanrichards")
  expect_equal(f$start[4], 1.3)
  expect_true(f$fix[4])
  expect_equal(unname(predict(f, newdata=data.frame(dbh=0))), 1.3)
  expect_equal(logLik(f), logLik(gcTMB(height ~ dbh, data=PineHeightDiameter, model="chapmanrichards", start=c(b4=1.3))))
  set.seed(5)
  dk <- data.frame(dbh=runif(150, 2, 60))
  dk$h <- gcTMB:::gcmeanR(mko, dk$dbh, 30, 14, 2, 1.3) + rnorm(150, 0, 1)
  w <- gcTMB(h ~ dbh, data=dk, model="korf")
  expect_equal(unname(predict(w, newdata=data.frame(dbh=0))), 1.3)
  expect_equal(w$start[4], 1.3)
  expect_equal(w$fix[3:4], c(FALSE, TRUE))
  # given in start, or estimated from the breast height
  expect_equal(unname(gcTMBstart(height ~ dbh, data=PineHeightDiameter, model="korf", start=c(b4=0))[4]), 0)
  expect_false(gcTMB(height ~ dbh, data=PineHeightDiameter, model="chapmanrichards", fix=c(b4=FALSE))$fix[4])
  # the threshold of the sigmoid curves and the start of beta
  for (m in c("logistic", "gompertz", "beta_capped")){
    fm <- gcTMB(height ~ DAS, data=FabaHeights25, model=m)
    expect_equal(fm$start[4], 0, label=m)
    expect_true(fm$fix[4], label=m)
  }
  # preecebaines estimates its rates b4 and b5
  fp <- gcTMB(height ~ age, data=subset(BerkeleyGrowth, sex == "female" & age >= 2), model="preecebaines")
  expect_false(any(fp$fix[4:5]))
})

test_that("the start b4 of beta is estimated, and the starting values are close to the estimates", {
  d <- simstart(mbeta, c(50, 28, 48, 3), sd=0.5)
  f <- gcTMB(y ~ day, data=d, model="beta_capped", fix=c(b4=FALSE))
  expect_equal(f$opt$convergence, 0)
  expect_equal(unname(coef(f)[1:4]), c(50, 28, 48, 3), tolerance=0.1)
  for (m in c("chapmanrichards", "beta_capped", "korf")){
    d <- simstart(gcTMB:::modelindex(m), switch(m, chapmanrichards=c(50, 0.06, 3, 0), beta_capped=c(50, 25, 45, 0), korf=c(50, 22, 2, 0)))
    sv <- gcTMBstart(y ~ day, data=d, model=m, start=c(b4=0))
    f0 <- gcTMB(y ~ day, data=d, model=m, start=c(b4=0))
    used <- 1:3
    expect_equal(unname(sv[used]), unname(f0$estimates[used]), tolerance=0.3, label=m)
  }
})

test_that("timeto of chapmanrichards, beta and korf: the inverse of the curve, with delta method standard errors", {
  for (p in list(list("chapmanrichards", c(50, 0.06, 3, 0), NULL), list("beta_capped", c(50, 25, 45, 0), NULL),
                 list("korf", c(50, 22, 2, 0), NULL), list("korf", c(50, 3.6, 1.6, 0), NULL))){
    mod <- gcTMB:::modelindex(p[[1]])
    d <- simstart(mod, p[[2]])
    f <- gcTMB(y ~ day, data=d, model=p[[1]], start=c(b4=0), fix=p[[3]])
    tt <- timeto(f, level=c(10, 50, 90))
    pv <- f$start
    pv[!f$fix] <- f$estimates[seq_len(sum(!f$fix))]
    gctimeR <- gcTMB:::gctimeR
    timefun <- function(par, q){
      pv[!f$fix] <- par
      gctimeR(mod, q, pv[2], pv[3], pv[4])
    }
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
    b1 <- gcTMB:::linkfun(gcTMB:::linkindex(f$link[["b1"]]), pv[1], inverse=TRUE)
    expect_equal(gcTMB:::gcmeanR(mod, coef(ta)[[1]], b1, pv[2], pv[3], pv[4]), 20, tolerance=1e-8)
    expect_true(all(is.finite(coef(maxrate(f))[1:2])))
  }
  # profile intervals of a time and of the inflection point
  for (m in c("chapmanrichards", "korf")){
    d <- simstart(gcTMB:::modelindex(m), if (m == "korf") c(50, 22, 2, 0) else c(50, 0.06, 3, 0))
    f <- gcTMB(y ~ day, data=d, model=m, start=c(b4=0))
    t50 <- coef(timeto(f, level=50))
    ci <- confint(timeto(f, level=50))
    expect_true(ci[1, 1] < t50 && t50 < ci[1, 2], label=m)
    mr <- maxrate(f)
    cim <- confint(mr, parm="tmax")
    expect_true(cim[1, 1] < coef(mr)[["tmax"]] && coef(mr)[["tmax"]] < cim[1, 2], label=m)
  }
})

test_that("b1 as the maximum growth rate of the chapmanrichards, beta and korf curves", {
  for (p in list(list("chapmanrichards", c(50, 0.06, 3, 0)), list("beta_capped", c(50, 25, 45, 0)), list("korf", c(50, 22, 2, 0)))){
    d <- simstart(gcTMB:::modelindex(p[[1]]), p[[2]])
    fa <- gcTMB(y ~ day, data=d, model=p[[1]], start=c(b4=0))
    fr <- gcTMB(y ~ day, data=d, model=p[[1]], start=c(b4=0), b1type="rate", link=c(b1="identity"))
    expect_equal(logLik(fr), logLik(fa), tolerance=1e-6, label=p[[1]])
    expect_equal(coef(fr)[[1]], coef(maxrate(fa))[["rate"]], tolerance=1e-5, label=p[[1]])
  }
})

# the beta curve not capped at its maximum: eqn (8) of Yin et al. (2003) beyond the end of growth b3,
# back to zero at 2 b3 - b2
mbu <- 12
betafull <- function(t, b1, b2, b3, b4){
  D <- b3 - b4
  d <- (b3 - b2)/D
  u <- (t - b4)/D
  ifelse(t > b4 & u < 1 + d, b1 * (1 + (1 - u)/d) * pmax(u, 0)^(1/d), 0)
}

test_that("the beta curve falls back to zero at 2 b3 - b2, and is the capped one up to b3", {
  t <- c(-5, 0, 5, 20, 25, 40, 45, 50, 60, 64.9, 65, 80)
  for (b4 in c(0, 3)){
    expect_equal(gcTMB:::gcmeanR(mbu, t, 50, 25, 45, b4), betafull(t, 50, 25, 45, b4))
    up <- t <= 45
    expect_equal(gcTMB:::gcmeanR(mbu, t[up], 50, 25, 45, b4), gcTMB:::gcmeanR(mbeta, t[up], 50, 25, 45, b4))
  }
  # the maximum b1 at b3, zero from 2 b3 - b2 = 65 on, and a negative growth rate after b3
  expect_equal(gcTMB:::gcmeanR(mbu, c(45, 65, 70), 50, 25, 45, 0), c(50, 0, 0))
  expect_true(all(diff(gcTMB:::gcmeanR(mbu, seq(45, 65, by=1), 50, 25, 45, 0)) < 0))
  # times to reach a level and the maximum growth rate are those of the rising part
  q <- c(0.1, 0.5, 0.9)
  expect_equal(gcTMB:::gctimeR(mbu, q, 25, 45, 0), gcTMB:::gctimeR(mbeta, q, 25, 45, 0))
  expect_equal(gcTMB:::maxrateR(mbu, c(50, 25, 45, 0, 1)), gcTMB:::maxrateR(mbeta, c(50, 25, 45, 0, 1)))
})

test_that("beta fits agree with nls, with timeto and maxrate from the TMB template", {
  set.seed(6)
  d <- data.frame(day=rep(seq(2, 64, by=3), 3))
  d$y <- betafull(d$day, 50, 25, 45, 0) + rnorm(nrow(d))
  f <- gcTMB(y ~ day, data=d, model="beta")
  expect_equal(f$opt$convergence, 0)
  n <- nls(y ~ betafull(day, b1, b2, b3, 0), data=d, start=list(b1=50, b2=25, b3=45))
  expect_equal(unname(bcoef(f)), unname(coef(n)), tolerance=1e-5)
  b <- bcoef(f)
  tt <- timeto(f, level=c(10, 50, 90))
  expect_equal(unname(coef(tt)), gcTMB:::gctimeR(mbu, c(0.1, 0.5, 0.9), b[2], b[3], 0), tolerance=1e-8)
  expect_equal(unname(coef(maxrate(f))), unname(gcTMB:::maxrateR(mbu, c(b, 0, 1))))
  # the capped curve does not fit the decline
  fc <- gcTMB(y ~ day, data=d, model="beta_capped")
  expect_true(AIC(f) < AIC(fc) - 10)
})
