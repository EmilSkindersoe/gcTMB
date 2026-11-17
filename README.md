gcTMB
================

The package gcTMB fits growth curves over time using the Template Model
Builder (TMB) framework, optionally with covariates on the parameters
and random effects.

## Installation

``` r
# You can install gcTMB from GitHub by
pak::pak("EmilSkindersoe/gcTMB")
```

## Documentation

The vignettes introduce the package with worked examples on faba bean
data:

1.  [Growth models and their parameters](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/growth-models.pdf): the
    growth models, the meaning of their parameters, confidence intervals,
    and simple treatment effects.
2.  [Fixed and random effects](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/fixed-random-effects.pdf):
    covariates and random effects, and likelihood ratio tests of their
    inclusion, for heights and counts.
3.  [Timing and maximum growth rates](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/timing-rates.pdf): `timeto()`
    and `maxrate()` with treatment and genotype effects, for normal and
    count data.
4.  [Autocorrelation and simulation](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/autocorrelation-simulation.pdf):
    autocorrelated residuals within individuals, and simulation studies.

[The overview](https://github.com/EmilSkindersoe/gcTMB/blob/main/vignettes/gcTMB.pdf) shows most of the features in one document.

## Growth models

The main function is `gcTMB()`. The models are of the form f(t) = b1
F(z(t)), with a function F increasing from 0 to 1 and a transformed time
z(t):

| `model` | f(t) | b2 | b3 |
|----|----|----|----|
| `"logistic"` | b1 / (1 + exp(-(t - b2)/b3)) | inflection point | scale |
| `"gompertz"` | b1 exp(-exp(-(t - b2)/b3)) | inflection point | scale |
| `"hill"` | b1 / (1 + (b2/t)^b3) | time of half the asymptote | slope on the log time axis |
| `"weibull"` | b1 (1 - exp(-(b3 - 1)/b3 (t/b2)^b3)) | inflection point | shape, \> 1 |
| `"normal"` | b1 Phi((t - b2)/b3) | inflection point | scale |
| `"gamma"` | b1 F(t), with F the CDF of the gamma distribution with shape b3 and scale b2/(b3 - 1) | inflection point | shape, \> 1 |
| `"linear"` | b1 t + b2, for comparison | intercept | not used |

b1 is the asymptote approached at late times, and the slope of the
linear model, which has no asymptote and is not truncated. The Hill,
Weibull and gamma curves are for positive times. The links of the
parameters can be chosen, e.g. `link = c(b1 = "log", b3 = "inverse")`,
with the identity (default), log, logit (b1) or inverse (b3 of the
logistic and Gompertz models, the rate 1/b3) links.

By default, the curve is truncated at a threshold time b4, (F(z(t)) -
F(z(b4)))/(1 - F(z(b4))), so that the response is zero up to b4. b4 is
fixed at 0, so that the curve starts at zero at t = 0, unless estimated
with `fix = c(b4 = FALSE)`; `truncated = FALSE` gives the curves without
truncation.

Each model has an asymmetry parameter b5 \> 0, fixed at 1 (the curves
above) unless estimated with `fix = c(b5 = FALSE)`: the logistic model
becomes the Richards curve b1 / (1 + b5 exp(-(t - b2)/b3))^(1/b5), the
Hill model b1 / (1 + (b2/t)<sup>b3)</sup>b5, the Weibull, normal and
gamma models a power b5 of F, and the Gompertz model b1 (1 - (1 -
exp(-exp(-(t - b2)/b3)))^b5), a power of its survival function.

Starting values come from a linearisation of the curve, Finv(y/b1) =
z(t), linear in t or log(t), profiled over the asymptote b1.

``` r
library(gcTMB)
set.seed(1)
d <- data.frame(day = rep(seq(0, 24, 2), each = 3))
d$weight <- 50 / (1 + exp(-(d$day - 10) / 2.5)) + rnorm(nrow(d), 0, 1.5)
fit <- gcTMB(weight ~ day, data = d, model = "logistic")
fit
#> Call:
#> gcTMB(form = weight ~ day, data = d, model = "logistic")
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>         50.099           9.770           2.647  
#> 
#> Residual standard deviation: 1.426
```

`timeto()` gives the times at which the curve reaches a percentage of
its asymptote, with standard errors from automatic differentiation and
profile likelihood intervals by `confint()`.

``` r
tt <- timeto(fit, level = c(10, 50, 90))
confint(tt)
#>         2.5 %    97.5 %
#> 10%  4.192012  4.902495
#> 50%  9.664311 10.137903
#> 90% 15.044105 16.338969
```
