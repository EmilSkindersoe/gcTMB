# Pairwise comparisons of a parameter between growth curves

Compares a parameter of the population curves (random effects at zero)
between every pair of rows of `newdata`, e.g. the treatment levels: the
difference of the parameter on the scale of its link, its standard error
from the covariance of the fixed effects, a Wald test of no difference,
and a confidence interval. With `b1type="rate"` in
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md), b1
is the maximum growth rate, and `pairwise(fit, newdata, "b1")` compares
the maximum growth rates.

## Usage

``` r
pairwise(x, newdata, parm = "b1", adjust = "tukey", level = 0.95, nsim = 1e+05)

# S3 method for class 'pairwise.gcTMB'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- newdata:

  Data frame. The values of the covariates in fform, one row for each
  curve to compare; the rows are named as in
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).

- parm:

  Character string. The parameter to compare, "b1" to "b5"

- adjust:

  Character string. Adjustment for the number of pairs: "tukey"
  (default) for the single-step adjustment of the p-values and
  simultaneous confidence intervals, a method of
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for adjusted
  p-values only, e.g. "holm", or "none"

- level:

  Numeric between 0 and 1. The level of the confidence intervals

- nsim:

  Integer. Number of simulations of the largest absolute test statistic
  for `adjust="tukey"`

- digits:

  Integer. Number of significant digits

- ...:

  Not used.

## Value

A data frame of class pairwise.gcTMB with a row for each pair: the
estimated difference, its standard error, the test statistic, the
p-value and the confidence limits, with the ratio and its limits for the
log link.

## Details

The tests use t quantiles with small-sample inference
([`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)),
and normal quantiles otherwise. By default (`adjust="tukey"`), the
p-values and the confidence intervals are adjusted by the single-step
method, the generalisation of Tukey's honestly significant differences
to any covariance of the estimates (Hothorn, Bretz and Westfall, 2008,
the default of the multcomp package): they use the distribution of the
largest absolute test statistic of all the pairs, given the correlations
of the differences, so that the intervals cover all the differences
simultaneously with probability `level`. This distribution is simulated,
with `nsim` draws, which gives the adjusted p-values to about 0.002 for
the default; the simulation has its own seed and leaves the random
number generator as it was. The other methods of
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html), e.g. "holm", adjust
the p-values only, and the confidence intervals are then not adjusted.
For a parameter with the log link, the differences are log ratios, and
the ratios with their confidence limits are given as well. Whether all
the curves have the same value of the parameter is tested by a
likelihood ratio test of the models with and without the covariates of
the parameter, see
[`anova.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/anova.gcTMB.md).

## References

Hothorn, T., Bretz, F. and Westfall, P. (2008). Simultaneous inference
in general parametric models. Biometrical Journal 50, 346-363.

## Examples

``` r
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
# b1 as the maximum growth rate, on the log scale: differences of log rates
fit <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25, b1type="rate", link=c(b1="log"))
pairwise(fit, newdata=data.frame(trt=levels(FabaHeights25$trt)), parm="b1")
#> Pairwise differences of b1 (maximum growth rate) on the scale of its log link:
#> 
#>                   Estimate Std.Error t value p.value   lower     upper ratio
#> trt=0.1 - trt=0.3  0.00614    0.0319   0.193  0.9997 -0.0817  0.093961 1.006
#> trt=0.1 - trt=0.5 -0.08075    0.0323  -2.502  0.0944 -0.1697  0.008194 0.922
#> trt=0.1 - trt=0.6 -0.08555    0.0327  -2.617  0.0719 -0.1756  0.004555 0.918
#> trt=0.1 - trt=0.8 -0.06441    0.0332  -1.940  0.3009 -0.1559  0.027119 0.938
#> trt=0.3 - trt=0.5 -0.08689    0.0312  -2.787  0.0462 -0.1728 -0.000941 0.917
#> trt=0.3 - trt=0.6 -0.09168    0.0316  -2.900  0.0340 -0.1788 -0.004529 0.912
#> trt=0.3 - trt=0.8 -0.07054    0.0322  -2.194  0.1881 -0.1592  0.018082 0.932
#> trt=0.5 - trt=0.6 -0.00479    0.0320  -0.150  0.9998 -0.0930  0.083414 0.995
#> trt=0.5 - trt=0.8  0.01634    0.0325   0.502  0.9876 -0.0733  0.106007 1.016
#> trt=0.6 - trt=0.8  0.02114    0.0329   0.642  0.9683 -0.0696  0.111918 1.021
#>                   ratio.lower ratio.upper
#> trt=0.1 - trt=0.3       0.922       1.099
#> trt=0.1 - trt=0.5       0.844       1.008
#> trt=0.1 - trt=0.6       0.839       1.005
#> trt=0.1 - trt=0.8       0.856       1.027
#> trt=0.3 - trt=0.5       0.841       0.999
#> trt=0.3 - trt=0.6       0.836       0.995
#> trt=0.3 - trt=0.8       0.853       1.018
#> trt=0.5 - trt=0.6       0.911       1.087
#> trt=0.5 - trt=0.8       0.929       1.112
#> trt=0.6 - trt=0.8       0.933       1.118
#> 
#> Wald t tests on 189 degrees of freedom; p-values and confidence limits adjusted by the single-step (Tukey) method, 95% confidence limits, and ratios.
```
