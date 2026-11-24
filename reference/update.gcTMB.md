# Update a gcTMB fit, or set its parameters

With `fixef`, `VarCorr`, `sigma` or `rho`, `update` returns the model
with some of its parameters set to new values, e.g. as the true model of
a simulation study (see
[`simulate.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md));
the model is not refitted. Without them, it refits the model with the
arguments in `...` changed, as
[`update`](https://rdrr.io/r/stats/update.html).

## Usage

``` r
# S3 method for class 'gcTMB'
update(object, ..., fixef = NULL, VarCorr = NULL, sigma = NULL, rho = NULL)
```

## Arguments

- object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- ...:

  Without fixef, VarCorr, sigma and rho: the arguments of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
  to change, as in [`update`](https://rdrr.io/r/stats/update.html).

- fixef:

  Named numeric vector or NULL. New values of coefficients of b1-b5,
  named as in
  [`fixef`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
  (e.g. "b1:(Intercept)"), or by the parameter (e.g. "b1") for a
  parameter with a single coefficient. Coefficients fixed in the fit,
  such as the asymmetry b5 (1 by default), can be set too, e.g.
  `fixef=c(b5=4)`; they stay fixed. NULL (default) keeps all the
  coefficients.

- VarCorr:

  List or NULL. New covariance matrices of the random effects, named by
  the terms as in
  [`VarCorr`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
  with the rows and columns in the same order; a single matrix for a
  model with one random effect term. The terms that are not given keep
  their covariance matrices, and NULL (default) keeps all of them.

- sigma:

  Numeric or NULL. The new residual standard deviation (gaussian, of
  log(y) for lognormal), or the dispersion parameter of the other
  families, as returned by
  [`sigma.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md).
  NULL (default) keeps it.

- rho:

  Numeric or NULL. The new autocorrelation rho of the residuals, as
  returned by
  [`corpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
  for a model fitted with a `correlation` structure: between 0 and 1 for
  [`corCAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md),
  between -1 and 1 for
  [`corAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md).
  NULL (default) keeps it.

## Value

A gcTMB object.

## Details

Every parameter of the model can be set: the coefficients of b1-b5
(`fixef`, including coefficients fixed in the fit), the covariance
matrices of the random effects (`VarCorr`), the residual standard
deviation or dispersion (`sigma`) and the autocorrelation (`rho`), i.e.
all the entries of `coef`. The parameters that are not given keep their
values, and the random effects are the conditional modes at the new
parameters. The new values are not estimates: the standard errors, the
covariance matrix of `vcov` and the confidence intervals are NA, and
[`converged`](https://emilskindersoe.github.io/gcTMB/reference/converged.md)
returns NA. The log-likelihood is that of the data at the new
parameters. To fit the data as well as possible with a chosen value of a
parameter, refit the model with that parameter fixed instead, e.g.
`correlation=corCAR1(0.8, form=~ time | ID, fixed=TRUE)` for the
autocorrelation, or `start=c(b5=2)` for the asymmetry, which is fixed by
default; the other parameters are then estimated given it, whereas
`update` keeps their values. For gaussian models without random effects,
the residual standard deviation is the value given, without the
small-sample correction of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

## Examples

``` r
fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
# larger asymptote, an asymmetric curve, and a larger variance of b1
vc <- VarCorr(fit)$Plant_ID
vc[1, 1] <- 2*vc[1, 1]
fit2 <- update(fit, fixef=c(b1=150, b5=4), VarCorr=list(Plant_ID=vc))
fixef(fit2)
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b5:(Intercept) 
#>      150.00000       30.86626       12.84146        4.00000 
VarCorr(fit2)
#> $Plant_ID
#>            b1        b2
#> b1 209.588060 -3.159913
#> b2  -3.159913  6.558684
#> attr(,"stddev")
#>        b1        b2 
#> 14.477156  2.560993 
#> attr(,"correlation")
#>             b1          b2
#> b1  1.00000000 -0.08522823
#> b2 -0.08522823  1.00000000
#> 
#> attr(,"sc")
#> [1] 2.074042
# refit with another growth model
fitg <- update(fit, model="gompertz")
# a weaker autocorrelation of the residuals, the other parameters unchanged
fitc <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25,
              correlation=corCAR1(form=~ DAS | Plant_ID))
corpars(update(fitc, rho=0.7))
#> rho 
#> 0.7 
```
