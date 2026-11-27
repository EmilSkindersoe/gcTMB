# Variance functions for heteroscedastic residuals

Variance functions for the `variance` argument of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md),
with the arguments and the parameterisation of `varIdent`, `varPower`,
`varConstPower`, `varExp` and `varComb` of the package nlme (whose
objects
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)
also accepts), and the polynomial `varPoly`, for the gaussian and
lognormal (on the log scale) families. The residual variance of an
observation is sigma^2 g^2, with the residual standard deviation sigma
of the model and the function g:

- `varIdent(form=~ 1 | s)`: g = delta_s, a ratio for each stratum s of
  the grouping variables, with delta = 1 for the reference stratum, so
  that each stratum has its own variance, e.g. a variance for each
  treatment (`~ 1 | Treatment`) or for each time point (`~ 1 | Time`); a
  numeric grouping variable is treated as a factor;

- `varPower(form=~ v)`: g = \|v\|^delta, variance sigma^2 \|v\|^(2
  delta);

- `varConstPower(form=~ v)`: g = c + \|v\|^delta, c \> 0, variance
  sigma^2 (c + \|v\|^delta)^2;

- `varExp(form=~ v)`: g = exp(t v), variance sigma^2 exp(2 t v);

- `varPoly(form=~ v, degree=k)`: g = 1 + d_1 (v - m) + d_2 (v - m)^2 +
  ... + d_k (v - m)^k, so that the standard deviation sigma g is a
  polynomial in v, e.g. one that is largest in the middle of the range
  of v (d_1 \> 0, d_2 \< 0, with its maximum at v = m - d_1/(2 d_2) for
  degree 2). The polynomial is centred at the mean m of v in the data
  (the mean response for `fitted(.)`, m = 0 for `rate(.)`), shown by
  print, so that sigma is the standard deviation at v = m, and the
  standard deviation at v = 0 can take any value; an uncentred
  polynomial, with sigma the standard deviation at v = 0, could not
  describe a standard deviation that would be 0 or below at v = 0 but is
  positive on the data, e.g. a straight line in the time from the first
  day. A polynomial can be 0 or negative, where the standard deviation
  is not defined: the likelihood is zero (the log-likelihood minus
  infinity) if the polynomial is 0 or below at any observation, so that
  the optimiser stays where it is positive on the data. It starts from d
  = 0, a constant standard deviation. Outside the data, e.g. in
  `predict` and `simulate` at new times, the polynomial can still be 0
  or below, which gives an error. A degree above 2 is rarely needed and
  tends to overfit;

- `varComb(...)`: the product of the multipliers g of two or more of the
  above, e.g. `varComb(varIdent(form=~ 1 | Treatment), varPower())` for
  a power of the fitted curve with a different factor for each
  treatment, or `varComb(varPower(), varExp())`, g = \|v\|^delta exp(t
  v), which with delta \> 0 and t \< 0 is largest at v = -delta/t, in
  the middle of growth. The variances are multiplied, not added: the
  standard deviations are sigma g_A g_B.

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

varExp(value = numeric(0), form = ~fitted(.), fixed = NULL)

varPoly(value = numeric(0), degree = 2, form = ~fitted(.), fixed = NULL)

varComb(...)
```

## Arguments

- value:

  Numeric. Starting values (or with `fixed`, see there): for `varIdent`
  the ratios of the standard deviations of the strata to the reference
  stratum, named by stratum (1 by default); for `varPower` the power and
  for `varExp` the exponent t (0 by default), named by level with a
  grouping factor; for `varPoly` the coefficients d_1, ..., d_k (0 by
  default), one for each power of v, or named `poly1`, ..., `polyk`
  (`poly1:<level>` for a single level of a grouping factor)

- form:

  One-sided formula: `~ 1 | s` for `varIdent`, with s one or more
  grouping variables (`~ 1 | A * B`); `~ v` or `~ v | g` for `varPower`,
  `varConstPower`, `varExp` and `varPoly`, with v a numeric variable,
  `fitted(.)` or `rate(.)`

- fixed:

  Values that are fixed and not estimated: for `varIdent` ratios named
  by stratum, for `varPower` the power and for `varExp` the exponent
  (named by level with a grouping factor), for `varConstPower` a list
  with `const` and/or `power`, e.g. `fixed=list(power=1)`, for `varPoly`
  coefficients named as in `value`, e.g. `fixed=c(poly2=0)`

- const, power:

  Numeric. For `varConstPower`, the starting values of the constant c
  (0.1 by default, as nlme) and the power delta (0 by default)

- degree:

  Integer. The degree k of the polynomial of `varPoly`, 2 by default

- ...:

  For `varComb`, two or more variance functions, whose standard
  deviation multipliers g are multiplied

## Value

An object of class `c("varIdent", "varFunc")`,
`c("varPower", "varFunc")`, `c("varConstPower", "varFunc")`,
`c("varExp", "varFunc")`, `c("varPoly", "varFunc")` or
`c("varComb", "varFunc")`, with the structure of the nlme objects (a
list of the components for `varComb`).

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
log ratios of `varIdent`), `var_par:power`, `var_par:const` (the log of
c), `var_par:expon` and `var_par:poly1`, `var_par:poly2`, ..., with
`:<level>` for a grouping factor, on the unconstrained scale of nlme;
[`varpars`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB-extract.md)
gives them on the natural scale. The parameters of the components of
`varComb` are prefixed with `A.`, `B.`, ..., as in nlme.

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
# \donttest{
# a standard deviation that is a quadratic polynomial of the fitted curve
fitp <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varPoly(degree=2))
varpars(fitp)
#>         poly1         poly2 
#>  8.109823e-03 -7.279771e-05 
# the product of a power and an exponential of the fitted curve
fitc <- gcTMB(Height ~ DAS, data=FabaHeights25, variance=varComb(varPower(), varExp()))
varpars(fitc)
#>      A.power      B.expon 
#>  0.833573792 -0.004632256 
# }
```
