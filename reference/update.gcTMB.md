# Set the parameters of a gcTMB fit

`update` returns the model with some of its parameters replaced by new
values, without fitting it again: the model keeps its structure (fform,
rform, model, family, correlation and variance), and

- `fixef` sets the coefficients of b1-b5 of `fform` (including
  coefficients fixed in the fit);

- `VarCorr` the covariance matrices of the random effects of `rform`;

- `sigma` the residual standard deviation, or the dispersion of the
  other families;

- `corpar` the parameters of the `correlation` structure (rho, or the
  range and nugget);

- `varpar` the parameters of the `variance` function.

The parameters that are not given keep their estimates, and the random
effects are the conditional modes at the new parameters. The result is
not a fit: it is moved away from the maximum of the likelihood, e.g. to
be the true model of a simulation study (see
[`simulate.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)).
The new values are not estimates: the standard errors, the covariance
matrix of `vcov` and the confidence intervals are NA, and
[`converged`](https://emilskindersoe.github.io/gcTMB/reference/converged.md)
returns NA. The log-likelihood is that of the data at the new
parameters.

## Usage

``` r
# S3 method for class 'gcTMB'
update(
  object,
  fixef = NULL,
  VarCorr = NULL,
  sigma = NULL,
  corpar = NULL,
  varpar = NULL,
  ...
)
```

## Arguments

- object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

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

- corpar:

  Numeric or NULL. New parameters of the correlation structure, as
  returned by
  [`corpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
  for a model fitted with a `correlation` structure: the autocorrelation
  rho, between 0 and 1 for
  [`corCAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  and between -1 and 1 for
  [`corAR1`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md);
  for
  [`corExp`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)
  and
  [`corGaus`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md),
  a named vector with the new range and/or nugget, e.g.
  `c(range=10, nugget=0.2)`. NULL (default) keeps them.

- varpar:

  Named numeric vector or NULL. New parameters of the variance function
  (see
  [`varFunc`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md))
  on their natural scale, named as returned by
  [`varpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
  e.g. `c(power=1)` or `c(const=2, power=0.5)` for `varConstPower`, or
  the ratio of a stratum for `varIdent`; the parameters not given keep
  their values. NULL (default) keeps them.

- ...:

  Not used; other arguments give an error that points to
  [`refit`](https://emilskindersoe.github.io/gcTMB/reference/refit.md).

## Value

A gcTMB object.

## Details

To fit the model again with other components, e.g. a new `fform`,
`rform`, `correlation`, `variance`, `model` or `family`, use
[`refit`](https://emilskindersoe.github.io/gcTMB/reference/refit.md);
the other arguments of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
give an error here. To set parameters of a model with a new component,
refit it first and then set the parameters of the new fit, e.g.
`update(refit(fit, fform=b1 ~ trt), fixef=...)`, as the coefficients of
the new fform are not known before the refit.

To fit the data as well as possible with a chosen value of a parameter,
refit with that parameter fixed instead, e.g.
`refit(fit, correlation=corCAR1(0.8, form=~ time | ID, fixed=TRUE))` for
the autocorrelation, or `refit(fit, start=c(b5=2))` for the asymmetry,
which is fixed by default; the other parameters are then estimated given
it, whereas setting it keeps their values. For gaussian models without
random effects, the residual standard deviation set by `sigma` is the
value given, without the small-sample correction of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

## See also

[`refit`](https://emilskindersoe.github.io/gcTMB/reference/refit.md) to
fit the model again with other components.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
# a larger asymptote, an asymmetric curve, and a larger variance of b1
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
# a weaker autocorrelation of the residuals, the other parameters unchanged
fitc <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25,
              correlation=corCAR1(form=~ DAS | Plant_ID))
corpars(update(fitc, corpar=0.7))
#> rho 
#> 0.7 
# a residual standard deviation proportional to the fitted curve
fitv <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varPower())
varpars(update(fitv, varpar=c(power=1)))
#> power 
#>     1 
```
