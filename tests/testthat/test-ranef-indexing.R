# the indexing of random effects over terms, levels and coefficients: the linear model f(t) = b1 t + b2
# with random effects is a linear mixed model, so gcTMB and lme4 (by maximum likelihood) give the same
# likelihood, variance components and predicted random effects, matched by level name

reindexdata <- function(){
  set.seed(42)
  ng <- 12
  nh <- 7
  d <- expand.grid(t=1:8, g=factor(paste0("g", 1:ng)))
  d$h <- factor(paste0("h", sample(nh, nrow(d), TRUE)))          # crossed with g
  d$a <- factor(paste0("a", (as.integer(d$g) - 1) %/% 3 + 1))   # g nested in a
  d$x <- rnorm(nrow(d))
  gi <- as.integer(d$g)
  ug <- cbind(rnorm(ng, 0, 2), rnorm(ng, 0, 1))
  d$y <- 10 + ug[gi, 1] + rnorm(nh, 0, 1.5)[as.integer(d$h)] + rnorm(4, 0, 2)[as.integer(d$a)] +
    (0.5 + rnorm(ng, 0, 0.8)[gi])*d$x + (2 + ug[gi, 2] + 0.5*ug[gi, 1])*d$t + rnorm(nrow(d))
  d
}

test_that("random effects of several terms, nested and crossed factors and slopes agree with lme4", {
  skip_if_not_installed("lme4")
  d <- reindexdata()
  cases <- list(list(r=b1 + b2 ~ 1 | g, f=NULL, l=y ~ t + (1 + t | g)),
                list(r=list(b1 ~ 1 | g, b2 ~ 1 | h), f=NULL, l=y ~ t + (0 + t | g) + (1 | h)),
                list(r=b2 ~ 1 | a/g, f=NULL, l=y ~ t + (1 | a/g)),
                list(r=b2 ~ x | g, f=b2 ~ x, l=y ~ t + x + (1 + x | g)),
                list(r=list(b1 + b2 ~ 1 | g, b2 ~ 1 | h), f=NULL, l=y ~ t + (1 + t | g) + (1 | h)))
  # columns of lme4's ranef for the random coefficients of gcTMB: b1 is the slope in t, b2 the intercept
  lcol <- c(b1="t", b2="(Intercept)", "b2:x"="x")
  for (k in seq_along(cases)){
    cs <- cases[[k]]
    # the nested case is singular: the variance of the outer factor is zero, as for lme4
    f <- suppressWarnings(gcTMB(y ~ t, fform=cs$f, rform=cs$r, data=d, model="linear"))
    m <- suppressMessages(lme4::lmer(cs$l, data=d, REML=FALSE, control=lme4::lmerControl(optimizer="bobyqa")))
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance=1e-6, label=paste("logLik, case", k))
    rg <- ranef(f)
    rl <- lme4::ranef(m)
    for (tn in names(rg)){
      L <- rl[[tn]]
      G <- rg[[tn]]
      expect_false(is.null(L), label=paste("term", tn))
      for (cn in colnames(G)) expect_equal(G[rownames(L), cn], L[[lcol[[cn]]]], tolerance=1e-4, label=paste("case", k, tn, cn))
    }
    # predictions for new data with the levels in another order, with and without the random effects
    nd <- d[sample(nrow(d), 40), ]
    nd$t <- nd$t + 0.5
    expect_equal(unname(predict(f, newdata=nd)), unname(predict(m, newdata=nd)), tolerance=1e-5, label=paste("predict, case", k))
    expect_equal(unname(predict(f, newdata=nd, random=FALSE)), unname(predict(m, newdata=nd, re.form=NA)), tolerance=1e-5)
    # simulations with the fitted random effects have the fitted means
    expect_equal(unname(attr(simulate(f, seed=1, random=FALSE), "mean")[[1]]), unname(fitted(m)), tolerance=1e-5)
    expect_equal(unname(attr(simulate(f, seed=1, random=FALSE, newdata=nd), "mean")[[1]]), unname(predict(m, newdata=nd)), tolerance=1e-5)
  }
  # times of fitted groups, with their predicted random effects
  f <- gcTMB(y ~ t, rform=b1 + b2 ~ 1 | g, data=d, model="linear")
  cf <- coef(lme4::lmer(y ~ t + (1 + t | g), data=d, REML=FALSE))$g
  tt <- timeto(f, level=30, type="absolute", newdata=data.frame(g=c("g7", "g2")))
  expect_equal(unname(coef(tt)), (30 - cf[c("g7", "g2"), "(Intercept)"])/cf[c("g7", "g2"), "t"], tolerance=1e-5)
})

test_that("simulate uses the fitted model, not the last parameters evaluated by TMB", {
  # a fit stopped short of the optimum; a later optimisation of its objective finds better parameters,
  # which TMB keeps as last.par.best
  fit <- suppressWarnings(gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25, control=list(iter.max=3)))
  m1 <- attr(simulate(fit, seed=1, random=FALSE), "mean")[[1]]
  invisible(nlminb(fit$opt$par, fit$obj$fn, fit$obj$gr))
  m2 <- attr(simulate(fit, seed=1, random=FALSE), "mean")[[1]]
  expect_equal(m2, m1)
  expect_equal(m2, unname(predict(fit)), tolerance=1e-8)
  expect_equal(simulate(fit, seed=2)[[1]], {
    fit2 <- suppressWarnings(gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25, control=list(iter.max=3)))
    simulate(fit2, seed=2)[[1]]
  })
  # control sdreport=FALSE: the random effects are the conditional modes at the estimates
  fn <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25, control=list(sdreport=FALSE))
  fs <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
  expect_equal(fn$sdrl$par.random, fs$sdrl$par.random, tolerance=1e-6)
})
