# Extract the components of a gcTMB fit

gcTMB has its own generics fixef, ranef and VarCorr and does not need
nlme or lme4. They have the same names as the generics of nlme
(re-exported by lme4), and work together with them: whichever package is
attached last, fixef, ranef and VarCorr work both for gcTMB fits and for
the fits of nlme and lme4 (the gcTMB generics pass other objects on to
nlme, if it is installed). The results for gcTMB fits have the format of
lme4: `ranef` a list of data frames, one for each term, also for a
single term (nlme returns a data frame for a single level of grouping),
and `VarCorr` a list of covariance matrices (nlme returns a character
table).

## Usage

``` r
fixef(object, ...)

ranef(object, ...)

VarCorr(x, sigma = 1, ...)

# Default S3 method
fixef(object, ...)

# Default S3 method
ranef(object, ...)

# Default S3 method
VarCorr(x, sigma = 1, ...)

# S3 method for class 'gcTMB'
fixef(object, ...)

# S3 method for class 'gcTMB'
ranef(object, ...)

# S3 method for class 'gcTMB'
VarCorr(x, sigma = 1, ...)

# S3 method for class 'gcTMB'
sigma(object, ...)

corpars(object)

varpars(object)
```

## Arguments

- object, x:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md);
  for the generics `fixef`, `ranef` and `VarCorr`, a fitted model of any
  class.

- ...:

  Not used.

- sigma:

  Not used; for compatibility with the generic of nlme.

## Value

`fixef` returns the coefficients of the parameters b1-b5, named as in
`coef`: the estimated coefficients and the coefficients fixed by the
user (e.g. an asymmetry b5 fixed at a value other than 1), without the
variance parameters and without the defaults b4 = 0 and b5 = 1.

`ranef` returns a list with a data frame for each random effect term,
named by its grouping factor: a row for each level and a column for each
random coefficient (e.g. b1, or b1:x for a random slope of b1 in x). The
random effects are the conditional modes, on the scale of the linear
predictors of the parameters (the scale of the link).

`VarCorr` returns a list with the covariance matrix of the random
effects of each term, named as in `ranef`, with the standard deviations
and the correlation matrix as attributes `stddev` and `correlation`; the
residual standard deviation is the attribute `sc`.

`sigma` returns the residual standard deviation of a gaussian model (of
log(y) for lognormal; the maximum likelihood estimate, while
[`nlme::gnls`](https://rdrr.io/pkg/nlme/man/gnls.html) reports it with
n - p degrees of freedom), or the dispersion parameter of the other
families as printed by `print.gcTMB` (precision, size or dispersion),
and NA for poisson and binomial models. The residual variance is
`sigma(object)^2`. For gaussian models without random effects and
small-sample inference (see
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)),
it is the unbiased estimate.

`corpars` returns the parameter of the autocorrelation of the residuals
within groups (see
[`corStruct`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)),
rho, the correlation of residuals one unit of time apart, and an empty
vector for a model with independent residuals; it is
`coef(fit$modelStruct$corStruct, unconstrained=FALSE)` of nlme.
`varpars` returns the parameters of the variance function of the
residuals, as a named vector; the current models have a constant
residual variance, so it returns an empty vector, as the extractor for a
variance function once it is added.

## Details

Extractors for the estimated parts of a model: the fixed effects, the
predicted random effects, the covariance matrices of the random effects,
the residual standard deviation, and the parameters of the residual
autocorrelation and of the variance function.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25)
fixef(fit)
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) 
#>      121.08844       30.86626       12.84146 
head(ranef(fit)$Plant_ID)
#>             b1         b2
#> 6    9.7005206  0.9458512
#> 7    0.3047415  1.5134969
#> 8  -10.0012561 -0.8431203
#> 9    9.4344827  2.9126002
#> 10   3.5648711 -0.9288945
#> 11   1.8398292  2.0643038
VarCorr(fit)
#> $Plant_ID
#>            b1        b2
#> b1 104.794030 -3.159913
#> b2  -3.159913  6.558684
#> attr(,"stddev")
#>        b1        b2 
#> 10.236896  2.560993 
#> attr(,"correlation")
#>            b1         b2
#> b1  1.0000000 -0.1205309
#> b2 -0.1205309  1.0000000
#> 
#> attr(,"sc")
#> [1] 2.074042
sigma(fit)^2
#> [1] 4.301651
corpars(fit)
#> named numeric(0)
varpars(fit)
#> named numeric(0)
```
