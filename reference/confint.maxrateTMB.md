# Confidence intervals for the maximum growth rate and its timing

Wald, profile likelihood or parametric bootstrap intervals for the
estimates of
[`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md).

## Usage

``` r
# S3 method for class 'maxrateTMB'
confint(
  object,
  parm = NULL,
  level = 0.95,
  method = c("wald", "profile", "boot"),
  nsim = 500,
  seed = NULL,
  boot = NULL,
  ...
)
```

## Arguments

- object:

  A maxrateTMB object, the result of
  [`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md).

- parm:

  Character or integer vector, or NULL. Estimates, names or positions in
  `coef(object)`; NULL (default) for all of them

- level:

  Numeric between 0 and 1. The confidence level

- method:

  Character string. "wald" for Wald intervals (default), "profile" for
  profile likelihood intervals, or "boot" for parametric bootstrap
  percentile intervals

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
  [`profile.maxrateTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.maxrateTMB.md)
  or
  [`bootstrap`](https://emilskindersoe.github.io/gcTMB/reference/bootstrap.md).

## Value

a matrix with the lower and upper confidence limits

## Details

Profile intervals are available for `tmax` and `rate`, see
[`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md);
for `fraction`, and for `tmax` of a curve at its threshold, Wald
intervals are returned with a message. Wald intervals use t quantiles
with small-sample inference
([`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)).
Bootstrap fits whose curve has no inflection point are left out of the
bootstrap intervals.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, data=FabaHeights25)
mr <- maxrate(fit)
confint(mr)
#>               2.5 %     97.5 %
#> tmax     29.3846302 32.0753316
#> rate      2.3743556  2.7456798
#> fraction  0.4386425  0.4693352
confint(mr, parm=c("tmax", "rate"), method="profile")
#> Computing profile confidence intervals ...
#>         2.5 %    97.5 %
#> tmax 29.43846 32.182525
#> rate  2.38040  2.752212
```
