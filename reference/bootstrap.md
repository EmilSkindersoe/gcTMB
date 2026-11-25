# Parametric bootstrap of a gcTMB model

Simulates responses from the fitted model, with new random effects (see
[`simulate.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)),
and fits the model again to each simulated response, with the iteration
limits of the fit (`control` in
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)),
starting at the estimates. The bootstrap sample can be used for
confidence intervals of the model parameters and of the times of
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md),
see
[`confint.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md),
without fitting the models again.

## Usage

``` r
bootstrap(object, nsim = 500, seed = NULL, progress = interactive())

# S3 method for class 'bootstrap.gcTMB'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'bootstrap.gcTMB'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- nsim:

  Integer. Number of simulations

- seed:

  Integer or NULL. Seed for the random number generator; NULL (default)
  sets no seed: the draws continue from the current state of the random
  number generator

- progress:

  Logical. Show a progress bar

- x:

  A bootstrap.gcTMB object, the result of `bootstrap`.

- digits:

  Integer. Number of significant digits in the printed output

- ...:

  Not used.

- parm:

  Character or integer vector, or NULL. Parameters, names or positions
  in `coef`; NULL (default) for all of them

- level:

  Numeric between 0 and 1. The confidence level

## Value

bootstrap returns a bootstrap.gcTMB object, a list with the estimates of
`coef(object)` for each simulation (`estimates`, a matrix with NA rows
for fits that did not converge), the fitted model (`fit`), and `nsim`
and `seed`. confint returns percentile intervals of the parameters, on
the scale of
[`confint.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.gcTMB.md).

## Examples

``` r
fit <- gcTMB(Height ~ DAS, data=FabaHeights25)
bt <- bootstrap(fit, nsim=50, seed=1, progress=FALSE)
bt
#> Parametric bootstrap: 50 of 50 fits converged
#> 
#>                Estimate    Mean     Bias Std. Error
#> b1:(Intercept)  120.785 120.939 0.153923     3.2508
#> b2:(Intercept)   30.730  30.760 0.029935     0.6897
#> b3:(Intercept)   12.881  12.927 0.046524     0.8880
#> log_sigma         2.077   2.081 0.003391     0.0522
confint(bt)
#>                     2.5 %     97.5 %
#> b1:(Intercept) 115.486761 126.788667
#> b2:(Intercept)  29.413291  32.173617
#> b3:(Intercept)  11.447898  14.613563
#> sigma            7.353607   8.738824
```
