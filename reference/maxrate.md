# Maximum growth rate and its timing

Estimates the time at which growth curves grow fastest, the inflection
point, the maximum growth rate there, and the fraction of the asymptote
reached at that time, with standard errors by the delta method.
Confidence intervals by `confint`: Wald, profile likelihood or
parametric bootstrap intervals.

## Usage

``` r
maxrate(x, curves = NULL, newdata = NULL, random = TRUE)
```

## Arguments

- x:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- curves:

  Numeric matrix, list of matrices, or NULL. Linear functions to
  identify the parameters of the growth curves, with 5 rows (b1-b5) and
  one column for each model coefficient, see
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).
  If curves and newdata are both NULL (default), the identity matrix,
  for models without covariates in fform.

- newdata:

  Data frame or NULL. The values of the covariates in fform, one row for
  each growth curve, as an alternative to curves; NULL (default) to use
  curves. With grouping variables of the random effects, the curves of
  fitted groups, see Details. The estimates are named by the covariate
  and group values (or by the row names of newdata, if set).

- random:

  Logical. If TRUE (default), the curves of rows of newdata with a level
  of a grouping variable in the fitted data include the predicted random
  effects of that group. If FALSE, all the curves are population curves,
  and the grouping variables in newdata are not used.

## Value

A maxrateTMB object: a list with the estimates (`coef`), their
covariance matrix (`vcov`), and what `confint` and `profile` need. The
estimates are named `tmax`, `rate` and `fraction`, followed by the name
of the curve for curves from `newdata` or named `curves`, e.g.
`rate:trt=0.1`; several unnamed `curves` are numbered.

## Details

The curves are the population curves, with all random effects at zero,
and the standard errors include the uncertainty of the fixed effects,
unless `newdata` names fitted groups.

If a row of newdata gives a level of a grouping factor that is in the
fitted data, e.g. `newdata=data.frame(Plant_ID="46_5")`, its curve is
the curve of that group, with its predicted random effects (the
conditional modes), as in
[`predict.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-methods.md)
and
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).
These are shrunken towards the population curve, more so for groups with
few observations. The standard errors then include the uncertainty of
the fixed effects and the prediction error of the random effects, from
the joint precision matrix of
[`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). A missing level
(NA), or a level that is not in the fitted data, gives the population
curve for that term, with a warning. Profile intervals are not available
for the curves of fitted groups; `confint` returns Wald intervals for
them, and bootstrap intervals with `method="boot"`.

For each curve, `maxrate` gives

- `tmax`: the time of the maximum growth rate, the inflection point of
  the curve;

- `rate`: the maximum growth rate, the slope of the curve at `tmax`, in
  units of the response per unit of time;

- `fraction`: the fraction of the asymptote b1 reached at `tmax`, how
  much of the growth has occurred when the growth is fastest.

For the logistic (Richards), hill and symmetric weibull and gompertz
curves these have closed forms; for the gompertz and weibull curves with
an asymmetry b5 different from 1 the inflection point is the root of an
equation, see the vignette. The hill and weibull curves have an
inflection point at a positive time only for b3 b5 \> 1; otherwise they
grow fastest at t = 0, with an unbounded rate for the hill curve, and
the results are NA. The linear model grows at the constant rate b1,
without an inflection point: its `tmax` and `fraction` are NA.

A truncated curve (the default of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md), at
b4 = 0) is the curve without truncation, f(t), shifted down to zero at
its threshold b4 and rescaled, b1 (f(t) - f(b4))/(b1 - f(b4)) beyond b4.
It has its maximum growth rate at the inflection point of the curve
without truncation if that is after the threshold, and otherwise just
after the threshold, where the slope jumps from zero: then `tmax` is b4
and `fraction` zero. For the logistic curve, truncated by default, the
inflection point is still at b2 (when b2 \> b4), but the maximum growth
rate is b1/(4 b3) divided by 1 - f(b4)/b1, larger than that of the curve
without truncation, and the fraction of the asymptote reached then is
below one half. The closed forms of the logistic curve without
truncation, `rate` = b1/(4 b3) and `fraction` = 1/2, hold only with
`truncated=FALSE`. The same holds for the curves of fitted groups, with
their own b1-b4.

The standard errors are by the delta method, with numerical derivatives
of the quantities with respect to the parameters of the curves. Profile
intervals are available for `tmax` (by b2, unless the curve is at its
threshold) and `rate` (by b1); for `fraction`, which depends on the
shape of the curve only, `confint` returns Wald intervals with a message
for `method="profile"`.

## See also

[`confint.maxrateTMB`](https://emilskindersoe.github.io/gcTMB/reference/confint.maxrateTMB.md)
for confidence intervals,
[`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)

## Examples

``` r
fit <- gcTMB(Height ~ DAS, data=FabaHeights25, model="logistic")
mr <- maxrate(fit)
mr
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>          Estimate Std.Error
#> tmax       30.730   0.68220
#> rate        2.560   0.09415
#> fraction    0.454   0.00778
confint(mr, method="wald")
#>               2.5 %     97.5 %
#> tmax     29.3846302 32.0753316
#> rate      2.3743556  2.7456798
#> fraction  0.4386425  0.4693352

# the Richards curve, with an estimated asymmetry
fitr <- gcTMB(Height ~ DAS, data=FabaHeights25, fix=c(b5=FALSE), link=c(b5="log"))
maxrate(fitr)
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>          Estimate Std.Error
#> tmax       30.666    3.0276
#> rate        2.560    0.0957
#> fraction    0.452    0.0854

# a curve for each treatment level
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fit2 <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25)
maxrate(fit2, newdata=data.frame(trt=levels(FabaHeights25$trt)))
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>                  Estimate Std.Error
#> tmax:trt=0.1       30.432   1.44817
#> rate:trt=0.1        2.441   0.10164
#> fraction:trt=0.1    0.452   0.00922
#> tmax:trt=0.3       29.080   1.42992
#> rate:trt=0.3        2.426   0.09967
#> fraction:trt=0.3    0.447   0.01029
#> tmax:trt=0.5       31.004   1.37572
#> rate:trt=0.5        2.647   0.10751
#> fraction:trt=0.5    0.454   0.00857
#> tmax:trt=0.6       31.507   1.39258
#> rate:trt=0.6        2.659   0.10806
#> fraction:trt=0.6    0.456   0.00821
#> tmax:trt=0.8       31.732   1.41979
#> rate:trt=0.8        2.604   0.10688
#> fraction:trt=0.8    0.457   0.00813

# the curves of two fitted plants, with their predicted random effects
fitre <- gcTMB(Height ~ DAS, rform=b1 ~ 1 | Plant_ID, data=FabaHeights25)
mrp <- maxrate(fitre, newdata=data.frame(Plant_ID=c(6, 7)))
mrp
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>                     Estimate Std.Error
#> tmax:Plant_ID=6       30.612   0.26353
#> rate:Plant_ID=6        2.715   0.05320
#> fraction:Plant_ID=6    0.453   0.00305
#> tmax:Plant_ID=7       30.612   0.26353
#> rate:Plant_ID=7        2.495   0.05099
#> fraction:Plant_ID=7    0.453   0.00305
confint(mrp)
#>                          2.5 %     97.5 %
#> tmax:Plant_ID=6     30.0957231 31.1287388
#> rate:Plant_ID=6      2.6103240  2.8188725
#> fraction:Plant_ID=6  0.4474963  0.4594599
#> tmax:Plant_ID=7     30.0957231 31.1287388
#> rate:Plant_ID=7      2.3949056  2.5947733
#> fraction:Plant_ID=7  0.4474963  0.4594599
```
