# Variance functions for heteroscedastic residuals

Variance functions for the `variance` argument of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
with the arguments and the parameterisation of `varIdent`, `varPower`
and `varConstPower` of the package nlme (whose objects
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
also accepts), for the gaussian and lognormal (on the log scale)
families. The residual variance of an observation is sigma^2 g^2, with
the residual standard deviation sigma of the model and the function g:

- `varIdent(form=~ 1 | s)`: g = delta_s, a ratio for each stratum s of
  the grouping variables, with delta = 1 for the reference stratum, so
  that each stratum has its own variance, e.g. a variance for each
  treatment (`~ 1 | Treatment`) or for each time point (`~ 1 | Time`); a
  numeric grouping variable is treated as a factor;

- `varPower(form=~ v)`: g = \|v\|^delta, variance sigma^2 \|v\|^(2
  delta);

- `varConstPower(form=~ v)`: g = c + \|v\|^delta, c \> 0, variance
  sigma^2 (c + \|v\|^delta)^2.

The covariate v of `varPower` and `varConstPower` is a numeric variable
of the data, e.g. `~ Time`, `fitted(.)` (the default), the fitted curve
f(t) of each observation (with its random effects), or `rate(.)`, the
growth rate f'(t) of the fitted curve at the time of the observation,
its slope, which is largest at the inflection point. `rate(.)` is a
gcTMB addition, not available in nlme; it is computed by a central
difference. With a grouping factor, `form=~ v | g`, each level of g has
its own power (and constant).

## Usage

``` r
varIdent(value = numeric(0), form = ~1, fixed = NULL)

varPower(value = numeric(0), form = ~fitted(.), fixed = NULL)

varConstPower(
  const = numeric(0),
  power = numeric(0),
  form = ~fitted(.),
  fixed = NULL
)
```

## Arguments

- value:

  Numeric. Starting values (or with `fixed`, see there): for `varIdent`
  the ratios of the standard deviations of the strata to the reference
  stratum, named by stratum (1 by default); for `varPower` the power (0
  by default), named by level with a grouping factor

- form:

  One-sided formula: `~ 1 | s` for `varIdent`, with s one or more
  grouping variables (`~ 1 | A * B`); `~ v` or `~ v | g` for `varPower`
  and `varConstPower`, with v a numeric variable, `fitted(.)` or
  `rate(.)`

- fixed:

  Values that are fixed and not estimated: for `varIdent` ratios named
  by stratum, for `varPower` the power (named by level with a grouping
  factor), for `varConstPower` a list with `const` and/or `power`, e.g.
  `fixed=list(power=1)`

- const, power:

  Numeric. For `varConstPower`, the starting values of the constant c
  (0.1 by default, as nlme) and the power delta (0 by default)

## Value

An object of class `c("varIdent", "varFunc")`,
`c("varPower", "varFunc")` or `c("varConstPower", "varFunc")`, with the
structure of the nlme objects.

## Details

The variance of `varPower` goes to zero where \|v\| does, e.g. the
fitted values of a curve that starts at zero or the growth rate before
and after growth, which makes those observations too influential and the
fit unstable; `varConstPower` keeps the variance at least sigma^2 c^2
there. A negative power gives an infinite variance at \|v\| = 0; with
`rate(.)`, this is the case for curves whose growth rate passes through
zero within the data, such as the beta curve after its maximum, which
the optimiser cannot handle, so `rate(.)` is meant for curves that keep
growing. Unlike nlme, which updates the fitted values between fits of
the variance parameters, gcTMB computes them inside the likelihood, so
the fit is the joint maximum likelihood estimate. The reference stratum
of `varIdent` is the first level of the grouping variables (sorted
numerically for a numeric variable), unless `value` and `fixed` name all
the others; nlme takes the first in the order of the data.

A variance function and a `correlation` structure can be combined: the
correlation is then that of the standardised residuals e/(sigma g). With
`varIdent(~ 1 | Time)`, `simulate` and prediction intervals need the
time points of the fit.

The estimated parameters are in `coef(fit)` as `var_par:<stratum>` (the
log ratios of `varIdent`), `var_par:power` and `var_par:const` (the log
of c), with `:<level>` for a grouping factor, on the unconstrained scale
of nlme;
[`varpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
gives them on the natural scale.

## References

Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and
S-PLUS. Springer, New York.

## See also

[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
[`varpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md),
[`corStruct`](https://emilskindersoe.github.io/gcTMB/reference/corStruct.md)

## Examples

``` r
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
# a residual variance for each treatment
fiti <- gcTMB(Height ~ DAS, fform=b1 ~ trt, data=FabaHeights25, variance=varIdent(form=~ 1 | trt))
varpars(fiti)
#>       0.1       0.3       0.5       0.6       0.8 
#> 1.0000000 0.9576577 0.6382966 0.2639973 0.3899443 
# residual standard deviation c + f'(t)^delta, largest where the curve grows fastest
fitr <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varConstPower(form=~ rate(.)))
varpars(fitr)
#>     const     power 
#>  2.834950 -2.224761 
```
