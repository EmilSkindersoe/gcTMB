# Profiling and benchmarking of gcTMB() against nlme
#
# The logistic model of sim_script_timed.R (the nlme project), with random effects on all three
# parameters and their full covariance matrix (pdLogChol), fitted with nlme and with gcTMB to the
# faba bean heights.

if (!requireNamespace("microbenchmark", quietly = TRUE)) install.packages("microbenchmark")
library(gcTMB)
library(nlme)
library(minpack.lm)
library(microbenchmark)

data_grouped <- groupedData(height ~ DAS | plant_ID, data = FabaHeights25)

# truncated logistic curve of the nlme project, with its gradient: A is the asymptote, tm the time of
# the maximum growth rate and k the scale; truncated at time 0 as gcTMB's default (b4 = 0)
lo <- function(x, A, tm, k, ...) {
  F_0 <- 1 / (1 + exp(tm / k))
  D   <- 1 - F_0
  F_x <- 1 / (1 + exp(-(x - tm) / k))
  N <- F_x - F_0
  value <- A * N / D
  d_A <- N / D
  dFx_dtm <- -F_x * (1 - F_x) / k
  dF0_dtm <- -F_0 * (1 - F_0) / k
  dFx_dk <- -(x - tm) * F_x * (1 - F_x) / k^2
  dF0_dk <-  tm * F_0 * (1 - F_0) / k^2
  d_tm <- A * ((dFx_dtm - dF0_dtm) * D + N * dF0_dtm) / D^2
  d_k  <- A * ((dFx_dk  - dF0_dk)  * D + N * dF0_dk)  / D^2
  attr(value, "gradient") <- cbind(A = d_A, tm = d_tm, k = d_k)
  value
}

# starting values from the logistic self-starter, as in the nlme project
start_logis <- getInitial(height ~ SSlogis(DAS, A, tm, k), data = data_grouped)
start_logis

# nlme, as fit_nlme_models() in the nlme project: nlsLM for starting values, then nlme with the
# control settings of the project
fit_nlme <- function() {
  init <- coef(nlsLM(height ~ lo(DAS, A, tm, k), data = data_grouped, start = start_logis))
  nlme(height ~ lo(DAS, A, tm, k),
       data    = data_grouped,
       fixed   = A + tm + k ~ 1,
       random  = pdLogChol(A + tm + k ~ 1),
       start   = init,
       control = nlmeControl(msMaxIter = 200, maxIter = 25, pnlsMaxIter = 20, minScale = 1e-8,
                             niterEM = 300, returnObject = FALSE, apvar = FALSE, opt = "nlminb"))
}

# gcTMB: the same model (b1 = A, b2 = tm, b3 = k, truncated at b4 = 0, symmetric b5 = 1), with
# correlated random effects on b1-b3, from gcTMB's own starting values (straight-line fits on the
# linearised curve), which give a better fit than the starting values of SSlogis
fit_gctmb <- function() {
  gcTMB(height ~ DAS, rform = b1 + b2 + b3 ~ 1 | plant_ID, data = FabaHeights25, model = "logistic")
}

# the two fits agree (the log-likelihoods differ a little: nlme approximates the likelihood by
# linearising the curve in the random effects, gcTMB by the Laplace approximation)
fn <- fit_nlme()
fg <- fit_gctmb()
rbind(nlme = c(fixef(fn), logLik = as.numeric(logLik(fn))),
      gcTMB = c(fixef(fg), logLik = as.numeric(logLik(fg))))

# benchmark, 3 repetitions each
bm <- microbenchmark(nlme = fit_nlme(), gcTMB = fit_gctmb(), times = 3)
print(bm, unit = "s")


## Profiling of gcTMB()

# 1. The steps of gcTMB(), timed one by one with the internal functions of the package, in the order
# gcTMB() calls them, over nrep repetitions
ns <- asNamespace("gcTMB")
rform <- b1 + b2 + b3 ~ 1 | plant_ID
nrep <- 5
steps <- c("data setup", "starting values", "MakeADFun", "optimisation", "sdreport")
times <- matrix(NA, nrep, length(steps), dimnames = list(NULL, steps))
evals <- matrix(NA, nrep, 2, dimnames = list(NULL, c("function", "gradient")))
for (r in seq_len(nrep)) {
  t0 <- proc.time()[["elapsed"]]
  tick <- function() { t1 <- proc.time()[["elapsed"]]; d <- t1 - t0; t0 <<- t1; d }
  dd <- ns$gcTMBdata(height ~ DAS, rform = rform, data = FabaHeights25, model = "logistic")
  fix <- ns$defaultfix(NULL, dd$npi)
  fix[dd$npi %in% ns$unusedpars(dd$mod, dd$trunc)] <- TRUE
  start <- ns$parvector(NULL, dd$npi, NA, "start")
  times[r, 1] <- tick()
  start <- ns$gcTMBstartvals(dd, start, fix, method = "lm")
  times[r, 2] <- tick()
  obj <- ns$gcTMBobj(dd, start, fix, rform = rform)
  times[r, 3] <- tick()
  ctrl <- list(iter.max = 1000, eval.max = 1000)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = ctrl)
  opt <- ns$retryvariances(obj, opt, dd$re$terms, -Inf, Inf, ctrl)
  times[r, 4] <- tick()
  sdr <- TMB::sdreport(obj, par.fixed = opt$par)
  times[r, 5] <- tick()
  evals[r, ] <- opt$evaluations
}
step_times <- rbind(mean = colMeans(times), share = colMeans(times) / sum(colMeans(times)))
round(step_times, 3)
evals[1, ]

# the time of one evaluation of the Laplace approximation (with the inner optimisation of the
# random effects) and of its gradient, at the estimates
microbenchmark(fn = obj$fn(opt$par), gr = obj$gr(opt$par), times = 20)

# the starting values: how long the variance starting values take within them
system.time(for (i in 1:5) ns$gcTMBstartvals(dd, ns$parvector(NULL, dd$npi, NA, "start"), fix, method = "lm"))

# 2. Rprof of whole fits, the functions with the largest total time (the compiled TMB code shows
# as .Call)
Rprof(prof_file <- tempfile(), interval = 0.005)
for (i in 1:5) fit_gctmb()
Rprof(NULL)
prof <- summaryRprof(prof_file)
head(prof$by.total, 25)
head(prof$by.self, 15)
