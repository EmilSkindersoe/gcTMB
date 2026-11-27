# Refit a gcTMB model with other components

`refit` fits a model again with some of its components changed, e.g. a
new `fform`, `rform`, `correlation`, `variance`, `model` or `family`:
the call of the original fit with the arguments in `...` changed and all
its other arguments reused (the formula, the data, the other structures,
`start`, `fix`, `control`, ...). The new fit is a fit like any other,
e.g. to compare with the original by
[`anova.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/anova.gcTMB.md)
or AIC. To change the parameters of a fit without fitting it again, e.g.
for the true model of a simulation, use
[`update.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/update.gcTMB.md).

## Usage

``` r
refit(object, ...)

# Default S3 method
refit(object, ...)

# S3 method for class 'gcTMB'
refit(object, ..., startfrom = c("data", "fit"))
```

## Arguments

- object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- ...:

  The arguments of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
  to change, named, e.g. `fform=b1 + b2 ~ trt`,
  `correlation=corCAR1(form=~ DAS | Plant_ID)` or `model="gompertz"`;
  NULL removes a component, e.g. `correlation=NULL`. An unnamed first
  argument changes the formula, as
  [`update.formula`](https://rdrr.io/r/stats/update.formula.html). For
  other objects, the arguments of the refit method of lme4.

- startfrom:

  Character string. The starting values: "data" (default) for those of
  [`gcTMBstart`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md),
  or "fit" for the estimates of `object` where they carry over, see
  Details.

## Value

A gcTMB object.

## Details

The starting values of the new fit are chosen by `startfrom`:

- `"data"` (default): from the data, as for a new call of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
  (see
  [`gcTMBstart`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBstart.md));

- `"fit"`: from the estimates of the fit being refitted, where they
  carry over, as a `start` list (see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)):
  the coefficients with the same name, growth model and link (e.g. the
  intercepts when covariates are added to `fform`, or all coefficients
  when only the correlation, the variance function or the random effects
  change); the variances and covariances of the random effects of the
  parameters on the same grouping factor (e.g. that of b1 when b2 is
  added to `b1 ~ 1 | Plant_ID`); the residual standard deviation or
  dispersion with the same family and variance function; and the
  parameters of the correlation structure and of the variance function
  when the refit keeps them (and they are not fixed). The rest starts
  from the data, e.g. the new coefficients of covariates, all
  coefficients of a new growth model, and the correlations of a random
  effect with a new one. With a new `rform`, the residual standard
  deviation of the fit may not suit the new random effects, but the
  optimiser moves it with them. A `start`, `correlation` or `variance`
  given in `...` is used as given. The call of the new fit shows the
  starting values used.

The data are found again by their name in the call, so a data frame
changed since the original fit gives a fit to the changed data.

The package lme4 has a `refit` function of its own, for a fit to a new
response. For objects that are not gcTMB fits, `refit` calls that of
lme4; if lme4 is attached after gcTMB, its `refit` masks this one, and
`gcTMB::refit` refits a gcTMB model.

## See also

[`update.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/update.gcTMB.md)
to set the parameters of a fit without refitting it.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
# another growth model, and covariates on b1 and b2, the other arguments reused
fitg <- refit(fit, model="gompertz")
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fitt <- refit(fit, fform=b1 + b2 ~ trt)
anova(fit, fitt)
#> Likelihood ratio tests of gcTMB models
#> 
#>      npar    AIC    BIC  logLik deviance     LR Df Pr(>LR)
#> fit     7 1041.3 1064.4 -513.67   1027.3                  
#> fitt   15 1046.3 1095.8 -508.17   1016.3 11.006  8  0.2014
# an autocorrelation structure with a random asymptote only
fitc <- refit(fit, rform=b1 ~ 1 | Plant_ID, correlation=corCAR1(form=~ DAS | Plant_ID))
# a refit that starts from the estimates of the fit it refits
fitc2 <- refit(fitc, fform=b1 + b2 ~ trt, startfrom="fit")
fitc2$call
#> gcTMB(form = Height ~ DAS, fform = b1 + b2 ~ trt, rform = b1 ~ 
#>     1 | Plant_ID, data = FabaHeights25, start = list(fixef = c(121.07040092474, 
#> NA, NA, NA, NA, 30.6368811016357, NA, NA, NA, NA, 12.9877198644004, 
#> 0, 1), VarCorr = list(Plant_ID = 106.828045118857), sigma = 4.30919732396573), 
#>     correlation = corCAR1(value = 0.97266882677351, form = ~DAS | 
#>         Plant_ID))
```
