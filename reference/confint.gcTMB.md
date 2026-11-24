# Confidence intervals for gcTMB models and times to reach a level of the asymptote

Profile likelihood intervals (default), Wald intervals, or parametric
bootstrap intervals for the coefficients and variance parameters of a
gcTMB model, or for the times estimated by
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).

## Usage

``` r
# S3 method for class 'gcTMB'
confint(
  object,
  parm = NULL,
  level = 0.95,
  method = c("profile", "wald", "boot"),
  nsim = 500,
  seed = NULL,
  boot = NULL,
  ...
)

# S3 method for class 'timetoTMB'
confint(
  object,
  parm = NULL,
  level = 0.95,
  method = c("profile", "wald", "boot"),
  nsim = 500,
  seed = NULL,
  boot = NULL,
  ...
)
```

## Arguments

- object:

  A gcTMB object
  ([`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md))
  or a timetoTMB object
  ([`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)).

- parm:

  Character or integer vector, or NULL. Parameters or times, names or
  positions in `coef(object)`; NULL (default) for all of them

- level:

  Numeric between 0 and 1. The confidence level

- method:

  Character string. "profile" for profile likelihood intervals, "wald"
  for Wald intervals, or "boot" for percentile intervals of a parametric
  bootstrap: responses are simulated from the fitted model (with new
  random effects, see
  [`simulate.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md))
  and the model is fitted again to each, starting at the estimates. Fits
  that do not converge are left out, with a message.

- nsim:

  Integer. Number of bootstrap simulations

- seed:

  Integer or NULL. Seed for the random number generator of the
  bootstrap; NULL (default) sets no seed: the draws continue from the
  current state of the random number generator

- boot:

  A bootstrap.gcTMB object or NULL. A bootstrap sample of the model from
  [`bootstrap`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md),
  used instead of new simulations; implies `method="boot"`. With NULL
  (default), method="boot" simulates nsim new bootstrap samples

- ...:

  Further arguments passed to
  [`profile.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md),
  [`profile.timetoTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.timetoTMB.md)
  or
  [`bootstrap`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md),
  e.g. `progress=FALSE`

## Value

a matrix with the lower and upper confidence limits

## Details

The counterpart of `intervals` of nlme, with profile likelihood and
parametric bootstrap intervals in addition to the Wald intervals; the
autocorrelation cor_par:AR1 or cor_par:CAR1 is reported as rho.

Standard deviations of the random effects, the residual standard
deviation and the precision are reported on their natural scale; Wald
intervals for them are computed on the log scale. The correlation of a
term with two random effects is reported as a correlation; the
parameters of larger correlation matrices on the unconstrained scale of
`coef`.

Profile and bootstrap intervals for the times are for the population
curves, and bootstrap intervals also for the curves of fitted groups
(`newdata` in
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)),
with the random effects of the groups drawn from their conditional
distribution at each bootstrap fit; for these, profile intervals are not
available and Wald intervals are returned with a message. Times computed
with `random=TRUE` for random effect terms without a fitted group, or
with `residual=TRUE`, include further variation in their standard
errors; for these, Wald intervals are returned with a message. Bootstrap
intervals need many fits of the model; with random effects this can take
several minutes. Computing profiles can take a while for models with
random effects; see
[`profile.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)
and
[`profile.timetoTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.timetoTMB.md),
whose results can be reused with `confint`.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, data=FabaHeights25)
confint(fit)
#> Computing profile confidence intervals ...
#>                     2.5 %    97.5 %
#> b1:(Intercept) 114.545439 129.48091
#> b2:(Intercept)  29.438468  32.18230
#> b3:(Intercept)  11.212413  15.18937
#> sigma            7.259687   8.83342
confint(fit, method="wald")
#>                    2.5 %     97.5 %
#> b1:(Intercept) 113.58225 127.987982
#> b2:(Intercept)  29.38463  32.075332
#> b3:(Intercept)  10.95416  14.807390
#> sigma            7.23722   8.804241
confint(timeto(fit, level=c(10, 50, 90)))
#> Computing profile confidence intervals ...
#>        2.5 %   97.5 %
#> 10% 10.02726 11.69745
#> 50% 31.32441 35.05932
#> 90% 55.61967 66.81695
# \donttest{
# parametric bootstrap, reused for the times
bt <- bootstrap(fit, nsim=200, seed=1, progress=FALSE)
confint(fit, boot=bt)
#>                     2.5 %     97.5 %
#> b1:(Intercept) 115.095303 128.532308
#> b2:(Intercept)  29.427031  32.151070
#> b3:(Intercept)  11.366012  14.939396
#> sigma            7.122528   8.655225
confint(timeto(fit, level=c(10, 50, 90)), boot=bt)
#>        2.5 %   97.5 %
#> 10% 10.19451 11.56292
#> 50% 31.43054 34.91577
#> 90% 55.90497 65.85907
# }
```
