# gcTMB Methods

gcTMB Methods

## Usage

``` r
# S3 method for class 'gcTMB'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'gcTMB'
summary(object, ...)

# S3 method for class 'summary.gcTMB'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'gcTMB'
logLik(object, ...)

# S3 method for class 'gcTMB'
nobs(object, ...)

# S3 method for class 'gcTMB'
df.residual(object, ...)

# S3 method for class 'gcTMB'
fitted(object, random = TRUE, ...)

# S3 method for class 'gcTMB'
residuals(
  object,
  type = c("pearson", "response", "quantile", "normalized"),
  random = TRUE,
  ...
)

# S3 method for class 'gcTMB'
predict(
  object,
  newdata = NULL,
  random = TRUE,
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  boot = NULL,
  nsim = 1000,
  seed = NULL,
  type = c("response", "parameters"),
  scale = c("original", "link"),
  ...
)

# S3 method for class 'gcTMB'
coef(object, ...)

# S3 method for class 'gcTMB'
vcov(object, ...)

# S3 method for class 'gcTMB'
simulate(object, nsim = 1, seed = NULL, random = TRUE, newdata = NULL, ...)
```

## Arguments

- x, object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- digits:

  Integer. Number of significant digits in the printed output

- ...:

  Not used.

- random:

  Logical. If TRUE, the random effects are included in the predictions,
  fitted values and residuals (conditional residuals). TRUE corresponds
  to the innermost level of `predict` and `fitted` of nlme, and FALSE to
  `level=0` (`re.form=NA` of lme4). For simulate, if TRUE, new random
  effects are drawn for each simulation, otherwise the fitted random
  effects are kept. For prediction intervals, random=TRUE gives the
  interval of a new observation in the group of each row of newdata, and
  random=FALSE the interval of a new observation in a new group,
  including the variance of the random effects. Rows of newdata with a
  level that is not in the fitted data are treated as new groups for
  that random effect term.

- type:

  Character string. For residuals, the type of residuals: "pearson"
  (default), the response residuals divided by the standard deviation of
  the fitted distribution and multiplied by the square root of the
  weights; "response", the observed minus the fitted values, as
  proportions for a binomial response; or "quantile", randomised
  quantile residuals (Dunn and Smyth, 1996), which are standard normal
  under the model. The types "response", "pearson" and "normalized" are
  those of `residuals` of nlme. "normalized" (gaussian and lognormal)
  are the standardised residuals decorrelated by the autocorrelation of
  the fit (see
  [`corStruct`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)):
  each standardised residual minus rho^dt times the previous one of its
  group, divided by sqrt(1 - rho^(2 dt)), independent standard normal
  under the model, as in nlme; without autocorrelation they are the
  pearson residuals; for the other families they are the pearson
  residuals, with a warning. For the gaussian family, the quantile
  residuals are multiplied by the square root of the weights, as
  precision weights with standard deviation sigma/sqrt(w); for the other
  families they do not use the weights, as a weighted log-likelihood is
  not the density of a distribution of the observation. The quantile
  residuals are random for the count families; use set.seed for
  reproducible values. For the lognormal family, the pearson and
  quantile residuals are the standardised residuals of log(y), (log(y) -
  log(f))/sigma, multiplied by the square root of the weights. For
  predict, the type of prediction: "response" (default), the mean
  response at the times of newdata, or "parameters", the parameters
  b1-b3 of the curve of each row of newdata (with its covariates in
  fform and, with random=TRUE, the random effects of its groups);
  newdata then needs no time. The standard errors of the parameters are
  those of the estimates, by the delta method, and their confidence
  intervals are computed on the scale of the link and back-transformed;
  with random=TRUE, they include the uncertainty of the predicted random
  effects of the fitted groups. Prediction intervals of the parameters
  add the variance of the random effects of the terms without a fitted
  level, the parameters of a new group, and are not available with boot.

- newdata:

  Data frame or NULL. The data for the predictions, with the time and
  the covariates of the model. If NULL, the original data is used. For
  simulate, the times, covariates and grouping variables of the
  observations to simulate, e.g. new times or new individuals; each
  level of a grouping variable is one curve, with random effects shared
  by its rows. With random=TRUE all levels get new random effects; with
  random=FALSE the levels in the fitted data keep their fitted random
  effects and the other levels get new ones. For binomial responses the
  number of trials is taken from the response in newdata (one for a
  binary response).

- se.fit:

  Logical. If TRUE, predict returns the standard errors of the predicted
  means, by the delta method; with random=TRUE they include the
  uncertainty of the predicted random effects.

- interval:

  Character string. The type of interval of predict: "none",
  "confidence" for confidence intervals of the mean, computed on the log
  scale for poisson and lognormal and on the logit scale for binomial
  and beta responses, or "prediction" for prediction intervals of a new
  observation with weight 1. Gaussian and lognormal models without
  random effects use t quantiles with n - p degrees of freedom and the
  unbiased residual variance, see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).
  For the gaussian and lognormal (on the log scale) families without
  boot, the prediction intervals use the normal (or t) distribution, and
  the variance of the random effects of new groups is added by the delta
  method. Otherwise they are quantiles of new observations simulated at
  nsim draws from the normal approximation of the fixed and random
  effects, with the variance parameters kept at their estimates, or at
  the bootstrap fits (see boot); with new random effects for new groups;
  for a binomial response as proportions, with the number of trials from
  the response in newdata (one for a binary response).

- level:

  Numeric between 0 and 1. Confidence level of the intervals of predict

- boot:

  A bootstrap.gcTMB object or NULL. A bootstrap sample of the model from
  [`bootstrap`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md).
  predict then gives the standard deviation of the bootstrap predictions
  as se.fit, percentile confidence intervals, and prediction intervals
  simulated at the bootstrap fits, including their variance parameters.
  The bootstrap responses have new random effects, so for random=TRUE
  the random effects of the fitted groups are drawn, for each bootstrap
  fit, from the Laplace approximation of their conditional distribution
  given the original data at the parameters of that fit.

- nsim:

  Integer. Number of response vectors to simulate; for predict, the
  number of parameter draws for simulated prediction intervals

- seed:

  Integer or NULL. An optional seed for the random number generator

- scale:

  Character string. For type="parameters", the scale of the parameters:
  "original" (default) or "link", the scale of the coefficients, see the
  argument link of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

## Value

Estimated model coefficients, variance-covariance matrix, or
predictions. predict returns a vector of predicted means, or, with
se.fit=TRUE or an interval, a data frame with the predicted means (fit),
their standard errors (se.fit), and the limits of the intervals (lwr,
upr). With type="parameters", it returns a matrix of the parameters used
by the model, a row for each row of newdata, or, with se.fit=TRUE or an
interval, a data frame with the row of newdata (row), the parameter, and
fit, se.fit, lwr and upr. coef and vcov include the variance parameters
(unlike `coef` of nlme and lme4, coef does not return the coefficients
of each group; see `predict(type="parameters")`): for each random effect
term the log standard deviations (log_sd:b1\|group, or
log_sd:b1:x\|group for a random slope of b1 in x) and, for correlated
random effects, the parameters of the correlation matrix
(corr_par:b2,b1\|group, unconstrained elements of the lower triangle of
its Cholesky factor, by row); and the log residual standard deviation
(log_sigma, gaussian; of log(y) for lognormal), log precision (log_phi,
beta), log size (log_theta, nbinom2), log dispersion (log_phi, nbinom1)
or log precision (log_phi, betabinomial); and for an estimated
autocorrelation the parameter cor_par:AR1 or cor_par:CAR1 on the
unconstrained scale of nlme (`coef` of an nlme corStruct), see
[`corStruct`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md).
summary returns a summary.gcTMB object with the coefficient table
(Estimate, Std. Error), the variance components, log-likelihood, AIC and
BIC, sample sizes, a Pearson dispersion statistic for poisson and
binomial models, and convergence information. logLik returns the
log-likelihood at the estimates, with the number of free coefficients
and variance parameters as degrees of freedom; it is the maximum
likelihood, comparable with nlme fits with `method="ML"`, not with the
REML fits of [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) and
[`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html). df.residual
returns the residual degrees of freedom n - p of gaussian and lognormal
models without random effects (small-sample inference, see
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)),
and NULL otherwise; as the `df` argument of
[`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html) for the
times of
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md),
it gives multivariate t instead of normal quantiles. simulate returns a
data frame with nsim columns of simulated responses at the estimated
parameters (or those set by
[`update.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/update.gcTMB.md)),
a row for each observation of the fit or row of newdata; for a binomial
response given as cbind(successes, failures), the number of successes.
Its attribute `ranef` has the random effects of each simulation, a list
named like the columns with an element in the format of
[`ranef`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
(for random=FALSE, the fitted random effects of the fitted levels), and
its attribute `mean` the true means of the simulated responses, a data
frame like the simulations: the curve f(t) of each observation with its
random effects (the median for lognormal; for binomial and beta-binomial
responses the proportion).

## Examples

``` r
fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
summary(fit)
#> Call:
#> gcTMB(form = height ~ DAS, rform = b1 + b2 ~ 1 | plant_ID, data = FabaHeights25)
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 200, groups: plant_ID 25
#> Log-likelihood: -513.7 (df = 7), AIC: 1041, BIC: 1064
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b1:(Intercept)   121.09      2.257
#> b2:(Intercept)    30.87      0.542
#> b3:(Intercept)    12.84      0.250
#> 
#> Random effects:
#>  Groups   Name Std.Dev. Corr 
#>  plant_ID b1   10.237        
#>           b2    2.561   -0.12
#> 
#> Residual standard deviation: 2.074
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 8.4e-05
coef(fit)
#>          b1:(Intercept)          b2:(Intercept)          b3:(Intercept) 
#>             121.0884446              30.8662604              12.8414641 
#>      log_sd:b1|plant_ID      log_sd:b2|plant_ID corr_par:b2,b1|plant_ID 
#>               2.3259984               0.9403950              -0.1214161 
#>               log_sigma 
#>               0.7294995 
logLik(fit)
#> 'log Lik.' -513.6748 (df=7)
# the population curve with confidence and prediction intervals
nd <- data.frame(DAS=c(20, 40, 60))
predict(fit, newdata=nd, random=FALSE, interval="confidence")
#>         fit   se.fit       lwr       upr
#> 1  28.69590 1.021166  26.69445  30.69735
#> 2  77.60741 1.772894  74.13260  81.08221
#> 3 108.71077 1.958312 104.87254 112.54899
predict(fit, newdata=nd, random=FALSE, interval="prediction")
#>         fit   se.fit      lwr       upr
#> 1  28.69590 1.021166 18.03090  39.36089
#> 2  77.60741 1.772894 59.58726  95.62755
#> 3 108.71077 1.958312 88.99044 128.43109
# the curve of plant 6, and its parameters
predict(fit, newdata=data.frame(DAS=c(20, 40, 60), plant_ID=6))
#> [1]  29.42109  81.76640 116.58430
predict(fit, newdata=data.frame(plant_ID=6), type="parameters")
#>           b1       b2       b3
#> [1,] 130.789 31.81211 12.84146
head(residuals(fit, type="quantile"))
#> [1]  0.7935475  1.1835400  0.7464733  0.3406020 -0.7506480  0.8659088
dim(simulate(fit, nsim=5, seed=1))
#> [1] 200   5
# three new plants at new times
nd <- expand.grid(DAS=seq(10, 80, 10), plant_ID=c("new1", "new2", "new3"))
sims <- simulate(fit, nsim=2, seed=1, newdata=nd)
head(sims)
#>       sim_1     sim_2
#> 1  15.84609  8.366161
#> 2  39.14342 24.477550
#> 3  66.27427 47.961726
#> 4  89.36434 69.243148
#> 5 110.53706 84.621353
#> 6 118.26703 95.285112
# the simulated random effects of the first simulation, and the true curves
attr(sims, "ranef")$sim_1
#> $plant_ID
#>             b1         b2
#> new1  6.282544 -4.2577268
#> new2 -1.906770 -0.7764702
#> new3  8.620931  1.8089401
#> 
head(attr(sims, "mean"))
#>       sim_1     sim_2
#> 1  14.83514  9.183935
#> 2  37.61210 24.600568
#> 3  65.08007 45.680227
#> 4  89.99773 67.660289
#> 5 107.40156 84.962581
#> 6 117.45848 95.810595
```
