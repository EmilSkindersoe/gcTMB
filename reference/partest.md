# Tests of the differences of a parameter between growth curves

Compares a parameter of the population curves (random effects at zero)
between the rows of `newdata`, e.g. the treatment levels: every row
against a reference row (`contrasts="dunnett"`, the default), or every
pair of rows (`contrasts="tukey"`). This is the function for the
hypothesis tests of the differences between curves; the summary of the
fit
([`summary.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md))
gives only the estimates and their standard errors. For each comparison
it gives the difference of the parameter on the scale of its link, its
standard error from the covariance of the fixed effects, a Wald test of
no difference, and a confidence interval. With `b1type="rate"` in
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md), b1
is the maximum growth rate, and `partest(fit, newdata, "b1")` compares
the maximum growth rates.

## Usage

``` r
partest(
  x,
  newdata,
  parm = "b1",
  contrasts = c("dunnett", "tukey"),
  ref = 1,
  adjust = "single-step",
  level = 0.95,
  nsim = 1e+05
)

# S3 method for class 'partest.gcTMB'
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

- contrasts:

  Character string. The comparisons: "dunnett" (default) for each row of
  newdata against the reference row `ref`, or "tukey" for all the pairs
  of rows

- ref:

  Integer or character string. The reference row of newdata for
  `contrasts="dunnett"`: its position, or its name as in the row names
  of the result; the first row by default. Not used for
  `contrasts="tukey"`

- adjust:

  Character string. Adjustment for the number of comparisons:
  "single-step" (default) for the single-step adjustment of the p-values
  and simultaneous confidence intervals (Tukey's or Dunnett's method, by
  `contrasts`), a method of
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for adjusted
  p-values only, e.g. "holm", or "none"

- level:

  Numeric between 0 and 1. The level of the confidence intervals

- nsim:

  Integer. Number of simulations of the largest absolute test statistic
  for `adjust="single-step"`

- digits:

  Integer. Number of significant digits

- ...:

  Not used.

## Value

A data frame of class partest.gcTMB with a row for each comparison,
named "a - b": the estimated difference a - b, its standard error, the
test statistic, the p-value and the confidence limits, with the ratio
and its limits for the log link.

## Details

The tests use t quantiles with small-sample inference
([`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)),
and normal quantiles otherwise. The contrasts and the adjustment for the
number of comparisons are separate choices, as in the multcomp package
(Hothorn, Bretz and Westfall, 2008): the contrasts say which differences
are tested, and the adjustment how their p-values and confidence
intervals account for the number of them.

By default (`adjust="single-step"`), the p-values and the confidence
intervals are adjusted by the single-step method, the default adjustment
of `summary.glht` in multcomp: they use the distribution of the largest
absolute test statistic of all the comparisons, given the correlations
of the differences, so that the intervals cover all the differences
simultaneously with probability `level`. For the comparisons with a
reference, this is Dunnett's test, and for all the pairs, Tukey's
honestly significant differences, both generalised to any covariance of
the estimates. Dunnett's test controls the familywise error rate of the
comparisons with the reference only, so its critical value is smaller
than Tukey's for the same rows. The distribution is simulated, with
`nsim` draws, which gives the adjusted p-values to about 0.002 for the
default; the simulation has its own seed and leaves the random number
generator as it was. The other methods of
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html), e.g. "holm", adjust
the p-values only, and the confidence intervals are then not adjusted;
"none" gives the unadjusted Wald tests.

For a parameter with the log link, the differences are log ratios, and
the ratios with their confidence limits are given as well. Whether all
the curves have the same value of the parameter is tested by a
likelihood ratio test of the models with and without the covariates of
the parameter, see
[`anova.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/anova.gcTMB.md).
With treatment contrasts in fform and the levels in newdata in the order
of the factor, the differences to the first level are the coefficients
of the parameter in the summary of the fit.

## References

Dunnett, C. W. (1955). A multiple comparison procedure for comparing
several treatments with a control. Journal of the American Statistical
Association 50, 1096-1121.

Hothorn, T., Bretz, F. and Westfall, P. (2008). Simultaneous inference
in general parametric models. Biometrical Journal 50, 346-363.

## Examples

``` r
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
# b1 as the maximum growth rate, on the log scale: differences of log rates
fit <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25, b1type="rate", link=c(b1="log"))
nd <- data.frame(trt=levels(FabaHeights25$trt))
# each treatment against the first level, by Dunnett's method
partest(fit, newdata=nd, parm="b1")
#> Differences to the reference of b1 (maximum growth rate) on the scale of its log link:
#> 
#>                   Estimate Std.Error t value p.value    lower  upper ratio
#> trt=0.3 - trt=0.1 -0.00614    0.0319  -0.193  0.9990 -0.08467 0.0724 0.994
#> trt=0.5 - trt=0.1  0.08075    0.0323   2.502  0.0454  0.00122 0.1603 1.084
#> trt=0.6 - trt=0.1  0.08555    0.0327   2.617  0.0336  0.00498 0.1661 1.089
#> trt=0.8 - trt=0.1  0.06441    0.0332   1.940  0.1655 -0.01744 0.1463 1.067
#>                   ratio.lower ratio.upper
#> trt=0.3 - trt=0.1       0.919        1.08
#> trt=0.5 - trt=0.1       1.001        1.17
#> trt=0.6 - trt=0.1       1.005        1.18
#> trt=0.8 - trt=0.1       0.983        1.16
#> 
#> Wald t tests on 189 degrees of freedom; p-values and confidence limits adjusted by the single-step (Dunnett) method, 95% confidence limits, and ratios.
# all the pairs, by Tukey's method
partest(fit, newdata=nd, parm="b1", contrasts="tukey")
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
