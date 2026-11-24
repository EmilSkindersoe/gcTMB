# Likelihood ratio tests for gcTMB models

Compares two or more models fitted to the same data by likelihood ratio
tests, e.g. a model with and without covariates in `fform`, with and
without random effects, or a poisson and a negative binomial model. The
test assumes nested models; the growth models (logistic, gompertz, hill,
weibull) are not nested, and are better compared by their AIC.

## Usage

``` r
# S3 method for class 'gcTMB'
anova(
  object,
  ...,
  test = c("LR", "F", "boot"),
  nsim = 1000,
  seed = NULL,
  progress = interactive()
)
```

## Arguments

- object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- ...:

  Further gcTMB objects, fitted to the same data.

- test:

  Character string. "LR" (default) for likelihood ratio tests, "F" for
  F-tests of gaussian or lognormal models without random effects, or
  "boot" for parametric bootstrap likelihood ratio tests. "Chisq", the
  name in earlier versions and in
  [`stats::anova`](https://rdrr.io/r/stats/anova.html), is accepted for
  "LR".

- nsim:

  Integer. Number of simulations for `test="boot"`

- seed:

  Integer or NULL. Seed for the random number generator, for
  `test="boot"`; NULL (default) sets no seed: the draws continue from
  the current state of the random number generator

- progress:

  Logical. Show a progress bar for `test="boot"`

## Value

a data frame of class anova with the number of parameters (npar), AIC,
BIC, log-likelihood and deviance (minus twice the log-likelihood) of
each model, and the likelihood ratio statistic (LR), its degrees of
freedom (Df) and p-value (Pr(\>LR)) for the comparison with the previous
model. For `test="F"`, the residual degrees of freedom (Res.Df) and
residual sum of squares (RSS) instead of the deviance, and the F
statistic (F) and its p-value (Pr(\>F)). For `test="boot"`, additionally
the bootstrap p-value (Pr(boot)), with the simulated likelihood ratio
statistics of each comparison in the attribute `"boot"`.

## Details

The models are ordered by their number of parameters, and each model is
compared with the previous one: the likelihood ratio statistic is twice
the difference of their log-likelihoods (the Laplace approximation of
the marginal log-likelihoods for models with random effects), with the
difference in the number of parameters as degrees of freedom. The test
assumes that the smaller model is nested in the larger one, which is not
checked. For parameters on the boundary of their range under the smaller
model, e.g. a standard deviation of random effects of zero, or an
infinite size of the negative binomial for the poisson model, the
chi-square distribution is conservative: the p-value is approximately
twice the correct value for a single standard deviation.

The tests are those of `anova` of nlme for models fitted with
`method="ML"`, in the layout of `anova` of lme4; as gcTMB fits by
maximum likelihood, models with different fixed effects can be compared
(REML fits of nlme cannot).

For gaussian and lognormal (on the log scale) models without random
effects, `test="F"` gives F-tests, as nls does for continuous responses:
the reduction of the (weighted) residual sum of squares per additional
parameter, divided by the residual mean square of the largest model,
with n - p residual degrees of freedom for p coefficients b1-b5. The
F-test is exact for linear models and accounts for the estimation of the
residual variance, which the likelihood ratio test does not, so it is
usually preferable for small experiments.

`test="boot"` gives parametric bootstrap likelihood ratio tests, as
`PBmodcomp` of the package pbkrtest: for each comparison, `nsim`
responses are simulated from the smaller model, with new random effects,
both models are fitted to each of them, and the p-value `Pr(boot)` is
the proportion of the simulated likelihood ratio statistics at least as
large as the observed one, (1 + \#{LR\* \>= LR}) / (1 + nsim). It does
not rely on the chi-square distribution, which is only approximate in
small samples, for parameters on the boundary (standard deviations of
random effects of zero, the poisson model within the negative binomial),
and for parameters that are not identified under the smaller model (e.g.
the correlation of random effects whose standard deviation is zero under
the smaller model). The chi-square p-value is shown for comparison. Each
model is fitted twice to each simulated response, starting at its
estimates and at the starting values of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
and the better fit is used; simulations without a finite log-likelihood
of both models are left out.

## Examples

``` r
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fit0 <- gcTMB(Height ~ DAS, data=FabaHeights25)
fit1 <- gcTMB(Height ~ DAS, fform=b1 ~ trt, data=FabaHeights25)
anova(fit0, fit1)
#> Likelihood ratio tests of gcTMB models
#> 
#>      npar    AIC    BIC  logLik deviance     LR Df  Pr(>LR)   
#> fit0    4 1406.5 1419.7 -699.23   1398.5                      
#> fit1    8 1396.2 1422.6 -690.10   1380.2 18.262  4 0.001097 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(fit0, fit1, test="F")
#> F-tests of gcTMB models
#> 
#>      npar    AIC    BIC  logLik Res.Df   RSS     F Df   Pr(>F)   
#> fit0    4 1406.5 1419.7 -699.23    197 12744                     
#> fit1    8 1396.2 1422.6 -690.10    193 11632 4.613  4 0.001404 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# random effects of the plants
fit2 <- gcTMB(Height ~ DAS, fform=b1 ~ trt, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25)
anova(fit1, fit2)
#> Likelihood ratio tests of gcTMB models
#> 
#>      npar    AIC    BIC  logLik deviance     LR Df   Pr(>LR)    
#> fit1    8 1396.2 1422.6 -690.10   1380.2                        
#> fit2    9 1134.2 1163.9 -558.09   1116.2 264.02  1 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
