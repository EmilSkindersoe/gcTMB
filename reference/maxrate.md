# Maximum growth rate and its timing

Estimates the time at which growth curves grow fastest, the inflection
point, the maximum growth rate there, and the fraction of the asymptote
reached at that time, with standard errors by the delta method.
Confidence intervals by `confint`: Wald, profile likelihood or
parametric bootstrap intervals.

## Usage

``` r
maxrate(x, curves = NULL, newdata = NULL)
```

## Arguments

- x:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- curves:

  Numeric matrix, list of matrices, or NULL. Linear functions to
  identify the parameters of the growth curves, with 5 rows (b1-b5) and
  one column for each model coefficient, see
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md);
  default is the identity matrix, for models without covariates in
  fform.

- newdata:

  Data frame or NULL. The values of the covariates in fform, one row for
  each growth curve, as an alternative to curves.

## Value

A maxrateTMB object: a list with the estimates (`coef`), their
covariance matrix (`vcov`), and what `confint` and `profile` need. The
estimates are named `tmax`, `rate` and `fraction`, followed by the name
of the curve for curves from `newdata` or named `curves`, e.g.
`rate:trt=0.1`; several unnamed `curves` are numbered.

## Details

The curves are the population curves, with all random effects at zero;
the standard errors include the uncertainty of the fixed effects. For
each curve, `maxrate` gives

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
and `fraction` zero.

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
fit <- gcTMB(height ~ DAS, data=FabaHeights25, model="logistic")
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
fitr <- gcTMB(height ~ DAS, data=FabaHeights25, fix=c(b5=FALSE), link=c(b5="log"))
maxrate(fitr)
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>          Estimate Std.Error
#> tmax       30.666    3.0276
#> rate        2.560    0.0957
#> fraction    0.452    0.0854

# a curve for each treatment level
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fit2 <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, data=FabaHeights25)
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
```
