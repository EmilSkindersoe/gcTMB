# Time to reach a level of the asymptote

Estimates the times at which growth curves reach a percentage of their
asymptote, or a given response, with standard errors by the delta method
and automatic differentiation.

## Usage

``` r
timeto(
  x,
  level = 50,
  curves = NULL,
  newdata = NULL,
  random = TRUE,
  residual = FALSE,
  type = c("relative", "absolute")
)
```

## Arguments

- x:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- level:

  Numeric vector. Levels in percent of the asymptote b1, between 0 and
  100 (default 50): the time at which the curve reaches p percent of b1,
  e.g. b2 + b3 log(p/(100 - p)) for the logistic curve, b2 - b3
  log(-log(p/100)) for gompertz, and b2 (p/(100 - p))^(1/b3) for hill.
  For the logistic, normal and hill curves without truncation, the time
  of 50 percent is b2. For a curve truncated at the threshold b4 (the
  default of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)),
  the levels are of the truncated curve, the time at which the curve f
  without truncation reaches f(b4) + p/100 (b1 - f(b4)). For the
  chapmanrichards and korf curves, the levels are of the growth above
  the intercept b4, and for preecebaines of the adult height b1. For
  `type="absolute"`, values of the response; the linear model, which has
  no asymptote, has only these.

- curves:

  Numeric matrix, list of matrices, or NULL. Linear functions to
  identify the parameters of the growth curve, with 5 rows (b1-b5) and
  one column for each model coefficient, or a list of such matrices. If
  curves and newdata are both NULL (default), the identity matrix, for
  models without covariates in fform. Each matrix is paired with an
  entry of level.

- newdata:

  Data frame or NULL. The values of the covariates in fform, one row for
  each growth curve, as an alternative to curves; the time is not
  needed. NULL (default) to use curves. With grouping variables of the
  random effects, the curves of fitted groups, see Details. The times
  are computed at each entry of level for each row, named by the
  covariate and group values (or by the row names of newdata, if set).

- random:

  Logical. If TRUE, the times of rows of newdata with fitted groups are
  conditional on their random effects, and for the other random effect
  terms, the variance-covariance matrix of the random effects is
  included in the standard error estimation. If FALSE, the times are for
  the population curves.

- residual:

  Logical. If TRUE, the residual variance is included in the standard
  error estimation (gaussian family only)

- type:

  Character string. "relative" (default) for levels relative to the
  asymptote, or "absolute" for the times at which the curves reach the
  responses in level, on the scale of the response (proportions for
  binomial, beta-binomial and beta). Responses outside the range of a
  curve, between zero and the asymptote, give NA, with a warning.

## Value

A timetoTMB object: a list with the estimated times (`coef`), their
covariance matrix (`vcov`), and what `confint` and `profile` need. The
times are named by the level, e.g. `10%` (relative) or `y=10`
(absolute), followed by the name of the curve for curves from `newdata`
or named `curves`, e.g. `10%:trt=0.1`; several unnamed `curves` are
numbered.

## Details

Without random effects, or with `random=FALSE`, the times are those of
the population curves, with all random effects at zero, and their
standard errors include the uncertainty of the fixed effects only. With
`random=TRUE`, the standard errors also include the variance of the
random effects, the variation of the time between groups.

If a row of newdata gives a level of a grouping factor that is in the
fitted data, the time of that row is for the curve of the group, with
its predicted random effects (the conditional modes), as in
[`predict.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md).
These are shrunken towards the population curve, more so for groups with
few observations. Their standard errors include the uncertainty of the
fixed effects and the prediction error of the random effects, from the
joint precision matrix of
[`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). Random effect
terms without a grouping variable in newdata, or with a level that is
not in the fitted data, are treated as above. Profile intervals are not
available for the times of fitted groups, see
[`confint.timetoTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md).

## See also

[`confint.timetoTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md)
for confidence intervals,
[`profile.timetoTMB`](https://emilskindersoe.github.io/gcTMB/reference/profile.timetoTMB.md)

## Examples

``` r
fit <- gcTMB(Height ~ DAS, data=FabaHeights25, model="logistic")
# days after sowing at which the plants reach 10%, 50% and 90% of their final height
tt <- timeto(fit, level=c(10, 50, 90))
tt
#> Times to reach the percentages of the asymptote:
#> 
#>     Estimate Std.Error
#> 10%     10.8     0.423
#> 50%     32.9     0.911
#> 90%     60.3     2.743
confint(tt)
#> Computing profile confidence intervals ...
#>        2.5 %   97.5 %
#> 10% 10.02726 11.69745
#> 50% 31.32441 35.05932
#> 90% 55.61967 66.81695
# the day at which they reach a height of 50
timeto(fit, level=50, type="absolute")
#> Times to reach the responses:
#> 
#>      Estimate Std.Error
#> y=50     28.8     0.336

# a curve for each treatment level
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fit2 <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25)
timeto(fit2, level=90, newdata=data.frame(trt=levels(FabaHeights25$trt)))
#> Times to reach the percentages of the asymptote:
#> 
#>             Estimate Std.Error
#> 90%:trt=0.1     60.2      2.90
#> 90%:trt=0.3     59.0      2.82
#> 90%:trt=0.5     60.7      2.92
#> 90%:trt=0.6     61.2      2.96
#> 90%:trt=0.8     61.4      2.98
```
