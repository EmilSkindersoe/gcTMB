# Profile likelihood of the maximum growth rate and its timing

Profiles the negative log-likelihood over the time of the maximum growth
rate (by b2) or the maximum growth rate (by b1) of a curve,
re-optimising the other coefficients.

## Usage

``` r
# S3 method for class 'maxrateTMB'
profile(
  fitted,
  parm = NULL,
  level = 0.99,
  maxit = 100,
  progress = interactive(),
  ...
)
```

## Arguments

- fitted:

  A maxrateTMB object, the result of
  [`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md).

- parm:

  Character or integer vector, or NULL. Estimates to profile, names or
  positions in `coef(fitted)`; tmax and rate only. NULL (default) for
  all of these.

- level:

  Numeric between 0 and 1. The highest confidence level the profiles are
  computed for.

- maxit:

  Integer. Maximum number of steps in each direction

- progress:

  Logical. Show a progress bar over the profiled estimates

- ...:

  Not used.

## Value

a profile.maxrateTMB object, which has the confint and plot methods of
[`profile.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)

## Examples

``` r
fit <- gcTMB(height ~ DAS, data=FabaHeights25)
pr <- profile(maxrate(fit), parm="rate")
confint(pr)
#>         2.5 %   97.5 %
#> rate 2.380395 2.752209
```
