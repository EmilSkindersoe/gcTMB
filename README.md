gcTMB
================

The package gcTMB fits growth curves over time using the Template Model Builder (TMB) framework, optionally with covariates on the parameters and random effects.

## Installation


``` r
# You can install gcTMB from GitHub by
pak::pak("EmilSkindersoe/gcTMB")
```

## Growth models

The main function is `gcTMB()`. The models are of the form f(t) = b1 F((t - b2)/b3), with the logistic (`model = "logistic"`) or Gompertz (`model = "gompertz"`) function F:

* b1 is the asymptote approached at late times,
* b2 is the time of the maximum growth rate (the inflection point),
* b3 is the scale of the curve on the time axis; with `link = c(b3 = "inverse")` its coefficients are the rate 1/b3.

The parameters b4 and b5 are reserved for future models. Starting values come from a linearisation of the curve, Finv(y/b1) = (t - b2)/b3, profiled over the asymptote b1.


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
#> Model:  logistic
#> Family: gaussian, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>         49.911           9.918           2.523  
#> 
#> Residual standard deviation: 1.363
```

`timeto()` gives the times at which the curve reaches a percentage of its asymptote, with standard errors from automatic differentiation and profile likelihood intervals by `confint()`.

``` r
tt <- timeto(fit, respLev = c(10, 50, 90))
confint(tt)
#>            2.5 %    97.5 %
#> 10:lf1  3.966569  4.771496
#> 50:lf2  9.697352 10.143677
#> 90:lf3 14.939651 16.016435
```
