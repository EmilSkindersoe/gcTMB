# Model selection criteria

Information criteria and goodness-of-fit criteria of one or more gcTMB
fits, for the choice between growth models fitted to the same data.

## Usage

``` r
criteria(..., which = NULL)
```

## Arguments

- ...:

  gcTMB objects, or a single list of gcTMB objects. The rows of the
  result are named by the names of the arguments or of the list, or
  numbered.

- which:

  Character vector or NULL. The criteria to return, in this order, e.g.
  `c("AIC", "R2", "BIC")`; NULL (default) returns all columns. A
  goodness-of-fit criterion without suffix, e.g. "R2", gives both its
  conditional and marginal version (R2.c and R2.m); "R2.m" gives only
  the marginal one.

## Value

A data frame with a row for each fit and the columns logLik, df, nobs,
the information criteria AIC, AICc, BIC, BICc, evBIC, HQIC, HQICc and
BC, and the goodness-of-fit criteria RMSE, ME, MAPE, RTSE, R2, R2adj,
CCC and CCCadj, each with suffix .c (conditional) and .m (marginal).

## Details

With L the maximised log-likelihood, k the number of estimated
parameters (free coefficients and variance parameters, the df of
`logLik`) and n the number of observations: AIC = -2L + 2k; AICc = AIC +
2k(k + 1)/(n - k - 1); BIC = -2L + log(n)k; BICc = BIC + k(k +
2)log(n)/(n - k - 2); evBIC = BIC + 2L/n; HQIC = -2L + 2k log(log(n));
HQICc = -2L + 2k log(log(n)) n/(n - k - 2); and the bridge criterion BC
= -2L + n^(2/3) (1 + 1/2 + ... + 1/k). All are smaller for better
models.

The goodness-of-fit criteria come in a conditional version (suffix .c),
with the curves of the groups including their random effects, and a
marginal version (suffix .m), with the population curve of the fixed
effects. The conditional versions use the residuals r = y - f(t) of the
observations; the marginal versions use the mean observation of each
time and population curve (the curves of the different values of the
covariates in fform), r = mean(y) - f(t), for each observation. RMSE is
the root mean squared residual, ME the mean residual (the bias), and
MAPE 100 times the mean absolute residual relative to the (mean)
observation. The root trapezoidal squared error RTSE approximates the
root mean integrated squared error between the observed and fitted
curves by the trapezoidal rule: the integral of r^2 over the times of
each curve divided by its time span, averaged over the curves, and the
square root. The curves are the levels of the random effect term with
the most levels; the marginal RTSE uses the population curves, and the
conditional RTSE is NA without random effects.

R2 is the coefficient of determination of Nakagawa and Schielzeth
(2013): with the variance of the marginal fitted values vf, the variance
of the conditional minus the marginal fitted values vr (zero without
random effects) and the variance of the conditional residuals ve, R2.m =
vf/(vf + vr + ve) and R2.c = (vf + vr)/(vf + vr + ve). CCC is the
concordance correlation coefficient of Vonesh et al. (1996) between the
observations and the conditional (CCC.c) or marginal (CCC.m) fitted
values. The adjusted versions R2adj and CCCadj are 1 - n/(n - p)(1 - R2)
with p the number of free coefficients. R2 and CCC are larger for better
models.

For binomial and beta-binomial responses the residuals are on the scale
of proportions.

## References

Nakagawa, S. and Schielzeth, H. (2013). A general and simple method for
obtaining R2 from generalized linear mixed-effects models. Methods in
Ecology and Evolution 4, 133-142.

Vonesh, E. F., Chinchilli, V. M. and Pu, K. (1996). Goodness-of-fit in
generalized nonlinear mixed-effects models. Biometrics 52, 572-587.

## Examples

``` r
fl <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
fg <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25, model="gompertz")
criteria(logistic=fl, gompertz=fg)
#>             logLik df nobs      AIC     AICc      BIC     BICc    evBIC
#> logistic -513.6748  7  200 1041.350 1041.933 1064.438 1066.185 1059.301
#> gompertz -530.0197  7  200 1074.039 1074.623 1097.128 1098.875 1091.827
#>              HQIC    HQICc       BC   RMSE.c    RMSE.m        ME.c        ME.m
#> logistic 1050.693 1051.793 1116.024 1.817786 0.5654044 -0.01070386 -0.02035558
#> gompertz 1083.383 1084.483 1148.714 2.055811 0.8278522  0.03554012  0.02100574
#>            MAPE.c   MAPE.m   RTSE.c    RTSE.m      R2.c      R2.m   R2adj.c
#> logistic 3.294648 0.651532 1.741021 0.6055063 0.9968727 0.9417571 0.9968250
#> gompertz 3.212344 1.422527 1.939751 0.8204698 0.9960090 0.9418983 0.9959482
#>            R2adj.m     CCC.c     CCC.m  CCCadj.c  CCCadj.m
#> logistic 0.9408701 0.9984232 0.9688412 0.9983992 0.9683667
#> gompertz 0.9410135 0.9979880 0.9686947 0.9979574 0.9682180
criteria(logistic=fl, gompertz=fg, which=c("AIC", "R2", "BIC"))
#>               AIC      R2.c      R2.m      BIC
#> logistic 1041.350 0.9968727 0.9417571 1064.438
#> gompertz 1074.039 0.9960090 0.9418983 1097.128
```
