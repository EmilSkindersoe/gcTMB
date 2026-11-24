# Likelihood profiles of times to reach a level of the asymptote

Profiles the (marginal) log-likelihood for the times estimated by
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).
At each value of a time, the model is reparameterised: the time, e.g.
b2 + b3 log(q/(1 - q)) for the logistic model, can be solved for b2, so
one coefficient of b2 is computed from the time and the other
coefficients, which are re-optimised. For truncated curves, the time
also depends on b2 through the threshold b4, and b2 is found
numerically.

## Usage

``` r
# S3 method for class 'timetoTMB'
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

  A timetoTMB object, the result of
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).

- parm:

  Character or integer vector, or NULL. Times to profile, names or
  positions in `coef(fitted)`; NULL (default) for all of them

- level:

  Numeric between 0 and 1. The highest confidence level the profiles are
  computed for.

- maxit:

  Integer. Maximum number of steps in each direction

- progress:

  Logical. Show a progress bar over the profiled times

- ...:

  Not used.

## Value

a profile.timetoTMB object, which has the confint and plot methods of
[`profile.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)

## Details

The profiles are for the times of the population curves, without the
random effects; the variation between groups (`random=TRUE` in
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md))
and the residual variation are not included. Gaussian models without
random effects use the F-based cutoff of
[`profile.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.gcTMB.md)
by default.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, data=FabaHeights25)
pr <- profile(timeto(fit, level=c(10, 90)), progress=FALSE)
confint(pr)
#>        2.5 %   97.5 %
#> 10% 10.02723 11.69743
#> 90% 55.61865 66.81623
plot(pr)

```
