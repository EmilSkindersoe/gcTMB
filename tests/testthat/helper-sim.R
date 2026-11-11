# Simulated data sets shared by the tests. All simulations use the gcTMB parameterisation,
# f(t) = b1 F((t - b2) / b3), via the mean function of the package.

models <- c("logistic", "gompertz")

# mean function for model name m and parameters b = c(b1, b2, b3)
gcmean <- function(m, t, b){
  gcTMB:::gcmeanR(match(m, models), t, b[1], b[2], b[3], 0, 0)
}

# inverse of F of each model
Finv <- list(logistic=function(q) qlogis(q), gompertz=function(q) -log(-log(q)))

# the same curves written out, for nls and nlme
formulas <- list(
  logistic=y ~ b1/(1 + exp(-(day - b2)/b3)),
  gompertz=y ~ b1*exp(-exp(-(day - b2)/b3))
)

# true parameters (asymptote, time of the maximum growth rate, scale) and observation times
truepars <- c(50, 10, 2.5)
days <- seq(0, 24, 2)

# gaussian data with nrep replicates per time
simgauss <- function(m, nrep=3, sd=1.5, seed=1, b=truepars, times=days){
  set.seed(seed)
  d <- data.frame(day=rep(times, each=nrep))
  d$y <- gcmean(m, d$day, b) + rnorm(nrow(d), 0, sd)
  d
}

# grouped data: ng groups (individuals), every time in every group, nrep replicates
simgroups <- function(times=days, ng=10, nrep=1, seed=1){
  set.seed(seed)
  d <- expand.grid(rep=seq_len(nrep), day=times, id=factor(seq_len(ng)))
  d$rep <- NULL
  d
}

# grouped gaussian data with random effects on b1 and b2
simre <- function(m="logistic", ng=12, sd1=5, sd2=1, sd=1, seed=3, b=truepars){
  d <- simgroups(ng=ng, seed=seed)
  u1 <- rnorm(ng, 0, sd1)
  u2 <- rnorm(ng, 0, sd2)
  d$y <- gcTMB:::gcmodel(match(m, models))$F((d$day - b[2] - u2[d$id])/b[3]) * (b[1] + u1[d$id]) + rnorm(nrow(d), 0, sd)
  d
}

# gcTMB coefficients without the variance parameters
bcoef <- function(fit){
  cf <- coef(fit)
  cf[grepl("^b[1-5]:", names(cf))]
}

# nls fit of model m to the faba bean heights of the package (FabaHeights25)
fabanls <- function(m, d=FabaHeights25){
  fo <- list(logistic=height ~ b1/(1 + exp(-(DAS - b2)/b3)),
             gompertz=height ~ b1*exp(-exp(-(DAS - b2)/b3)))[[m]]
  nls(fo, data=d, start=list(b1=120, b2=30, b3=15))
}

# nls fit of model m at the true parameters
nlsfit <- function(m, d, start=truepars){
  nls(formulas[[m]], data=d, start=list(b1=start[1], b2=start[2], b3=start[3]))
}
