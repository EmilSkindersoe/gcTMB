# Autocorrelation of the residuals within groups

Correlation structures for the residuals of the observations of a curve
over time (or in space), for the `correlation` argument of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
with the same arguments and the same parameterisation as `corAR1`,
`corCAR1`, `corExp` and `corGaus` of the package nlme (whose objects
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
also accepts).

## Usage

``` r
corAR1(value = 0, form = ~1, fixed = FALSE)

corCAR1(value = 0.2, form = ~1, fixed = FALSE)

corExp(
  value = numeric(0),
  form = ~1,
  nugget = FALSE,
  metric = c("euclidean", "maximum", "manhattan"),
  fixed = FALSE
)

corGaus(
  value = numeric(0),
  form = ~1,
  nugget = FALSE,
  metric = c("euclidean", "maximum", "manhattan"),
  fixed = FALSE
)
```

## Arguments

- value:

  Numeric. The starting value of the parameters, or their values with
  `fixed=TRUE`: rho, between -1 and 1 (default 0) for `corAR1` and
  between 0 and 1 (default 0.2) for `corCAR1`; for `corExp` and
  `corGaus` the range, or with a nugget the range and the nugget
  (default: see Details)

- form:

  One-sided formula. The time covariate (or the position variables,
  `~ x + y | group`, for `corExp` and `corGaus`) and the grouping
  factor, `~ time | group`. Without a grouping factor, the residuals of
  all observations form a single series when the model has no random
  effects, and are grouped by the innermost grouping factor of the
  random effects (the one with most levels) otherwise, as in
  [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html). Without a time
  covariate (`~ 1` or `~ 1 | group`), the times are the positions of the
  observations within their group.

- fixed:

  Logical. If TRUE, the parameters are fixed at `value` and not
  estimated

- nugget:

  Logical. If TRUE, the spatial structure has a nugget

- metric:

  Character. The distance between positions of several variables:
  "euclidean" (default), "maximum" or "manhattan", as in
  [`dist`](https://rdrr.io/r/stats/dist.html); all give \|t - s\| for a
  single time covariate

## Value

An object of class `c("corAR1", "corStruct")` or
`c("corCAR1", "corStruct")`: `value` on the unconstrained scale of nlme,
log((1 + rho)/(1 - rho)) for AR(1) and log(rho/(1 - rho)) for CAR(1),
with the attributes `formula` and `fixed`; or of class
`c("corExp", "corSpatial", "corStruct")` or
`c("corGaus", "corSpatial", "corStruct")`: `value` on the natural scale,
as in nlme, with the attributes `formula`, `nugget`, `metric` and
`fixed`.

## Details

The residuals of different groups are independent. Within a group, the
residuals e are stationary, with the residual variance sigma^2 of the
model (sigma^2/w for weights w). For `corAR1` and `corCAR1`, the
correlation is \$\$cor(e_s, e_t) = rho^{\|t - s\|},\$\$ with rho the
correlation of residuals one unit of time apart (Zuur et al., 2009,
Section 6.1; Pinheiro and Bates, 2000, Section 5.3):

- `corAR1`: the autoregressive process of order 1 at integer times
  (positions), with -1 \< rho \< 1; gaps in the times are allowed, e.g.
  measurements every second day with `form = ~ day | ID` and a day
  without measurement, with the correlation rho^2 across a gap of two
  days. Without a time covariate (`form = ~ 1 | ID`), the time is the
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

The spatial structures `corExp` and `corGaus` give the correlation as a
function of the distance d between the positions of two observations of
a group, with a range r \> 0 and, with `nugget=TRUE`, a nugget 0 \< n \<
1 (n = 0 otherwise):

- `corExp`: \\cor = (1 - n) \exp(-d/r)\\;

- `corGaus`: \\cor = (1 - n) \exp(-(d/r)^2)\\.

The position is the time covariate (`~ time | group`), with d the time
between the observations, or several covariates (`~ x + y | group`),
with d the distance of `metric`. Without a nugget, `corExp` in time is
`corCAR1` with rho = exp(-1/r). The nugget is the part of the residual
variance that is independent between observations, such as measurement
error, besides the correlated deviations of the process (Diggle, 1988);
with it, observations at the same position are allowed. The residuals of
a group are multivariate normal with this correlation matrix, so the
time of a fit grows with the cube of the number of observations per
group. As in nlme, the default starting value of the range is 90\\

The estimated parameters are in `coef(fit)` as `cor_par:AR1` or
`cor_par:CAR1`, or `cor_par:range` and `cor_par:nugget`, on the
unconstrained scale of nlme (`coef` of an nlme `corStruct`): the logit
link log(rho/(1 - rho)) for CAR(1), and log((1 + rho)/(1 - rho)), the
logit of (1 + rho)/2, for AR(1), whose rho lies between -1 and 1; log(r)
for the range and log(n/(1 - n)) for the nugget. Wald intervals are
computed on this scale and transformed back.
[`corpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
gives rho (or the range and nugget), and
`confint(fit, parm="cor_par:CAR1")` an interval for rho. If nlme is
attached after gcTMB, its constructors mask these; gcTMB accepts their
objects, and gcTMB's correlation structures can be passed to
[`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html),
[`nlme::gnls`](https://rdrr.io/pkg/nlme/man/gnls.html),
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) and
[`nlme::nlme`](https://rdrr.io/pkg/nlme/man/nlme.html).

## References

Diggle, P.J. (1988). An approach to the analysis of repeated
measurements. Biometrics 44, 959-971.

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
fit <- gcTMB(Height ~ DAS, fform=b1 + b2 ~ trt, rform=b1 + b2 ~ 1 | Plant_ID, data=FabaHeights25,
             correlation=corCAR1(form=~ DAS | Plant_ID))
corpars(fit)
#>       rho 
#> 0.9593852 
# spatial correlation between neighbouring trees of a plot, with a nugget (independent
# measurement error), on five of the plots
pines <- droplevels(subset(PineHeightDiameter, Plot %in% 1:5))
fitn <- gcTMB(Height ~ DBH, rform=b1 ~ 1 | Plot, data=pines, model="chapmanrichards",
              correlation=corExp(form=~ X + Y | Plot, nugget=TRUE))
corpars(fitn)
#>     range    nugget 
#> 8.4016652 0.2542902 
```
