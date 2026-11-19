# Autocorrelation of the residuals within groups

Correlation structures for the residuals of the observations of a curve
over time, for the `correlation` argument of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
with the same arguments and the same parameterisation as `corAR1` and
`corCAR1` of the package nlme (whose objects
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
also accepts).

## Usage

``` r
corAR1(value = 0, form = ~1, fixed = FALSE)

corCAR1(value = 0.2, form = ~1, fixed = FALSE)
```

## Arguments

- value:

  Numeric. The starting value of rho, or its value with `fixed=TRUE`:
  between -1 and 1 (default 0) for `corAR1`, between 0 and 1 (default
  0.2) for `corCAR1`

- form:

  One-sided formula. The time covariate and the grouping factor,
  `~ time | group`. Without a grouping factor, the residuals of all
  observations form a single series when the model has no random
  effects, and are grouped by the innermost grouping factor of the
  random effects (the one with most levels) otherwise, as in
  [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html). Without a time
  covariate (`~ 1` or `~ 1 | group`), the times are the positions of the
  observations within their group.

- fixed:

  Logical. If TRUE, rho is fixed at `value` and not estimated

## Value

An object of class `c("corAR1", "corStruct")` or
`c("corCAR1", "corStruct")`: `value` on the unconstrained scale of nlme,
log((1 + rho)/(1 - rho)) for AR(1) and log(rho/(1 - rho)) for CAR(1),
with the attributes `formula` and `fixed`.

## Details

The residuals of different groups are independent. Within a group, the
residuals e are stationary, with the residual variance sigma^2 of the
model (sigma^2/w for weights w) and the correlation \$\$cor(e_s, e_t) =
rho^{\|t - s\|},\$\$ with rho the correlation of residuals one unit of
time apart (Zuur et al., 2009, Section 6.1; Pinheiro and Bates, 2000,
Section 5.3):

- `corAR1`: the autoregressive process of order 1 at integer times
  (positions), with -1 \< rho \< 1; gaps in the times are allowed, e.g.
  measurements every second day with `form = ~ day | plant` and a day
  without measurement, with the correlation rho^2 across a gap of two
  days. Without a time covariate (`form = ~ 1 | plant`), the time is the
  position of the observation within its group, in the order of the
  data;

- `corCAR1`: the continuous-time autoregressive process of order 1, for
  times that need not be integers and unequally spaced observations,
  with 0 \< rho \< 1.

The times must be unique within a group. As the process is Markov, the
likelihood is the product of the density of the first residual of each
group and the conditional densities of each residual given the previous
one, normal with mean rho^dt e_prev and variance sigma^2 (1 - rho^(2
dt)) for the time dt between them, which is exact and linear in the
number of observations.

The estimated parameter is in `coef(fit)` as `cor_par:AR1` or
`cor_par:CAR1`, on the unconstrained scale of nlme (`coef` of an nlme
`corStruct`): the logit link log(rho/(1 - rho)) for CAR(1), and log((1 +
rho)/(1 - rho)), the logit of (1 + rho)/2, for AR(1), whose rho lies
between -1 and 1. Wald intervals are computed on this scale and
transformed back.
[`corpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
gives rho, and `confint(fit, parm="cor_par:CAR1")` an interval for rho.
If nlme is attached after gcTMB, its constructors mask these; gcTMB
accepts their objects, and gcTMB's correlation structures can be passed
to [`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html),
[`nlme::gnls`](https://rdrr.io/pkg/nlme/man/gnls.html),
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) and
[`nlme::nlme`](https://rdrr.io/pkg/nlme/man/nlme.html).

## References

Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and
S-PLUS. Springer, New York.

Zuur, A.F., Ieno, E.N., Walker, N.J., Saveliev, A.A. and Smith, G.M.
(2009). Mixed Effects Models and Extensions in Ecology with R. Springer,
New York.

## See also

[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
[`corpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
for the estimate of rho, and `residuals(fit, type="normalized")` for
residuals that are independent under the model.

## Examples

``` r
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
fit <- gcTMB(height ~ DAS, fform=b1 + b2 ~ trt, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25,
             correlation=corCAR1(form=~ DAS | plant_ID))
corpars(fit)
#>       rho 
#> 0.9593852 
```
