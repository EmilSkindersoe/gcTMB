# Starting values for a gcTMB model

Computes starting values for the fixed-effect parameters of a growth
curve, which can be inspected, modified, and passed to the `start`
argument of
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

## Usage

``` r
gcTMBstart(
  form,
  fform = NULL,
  data,
  start = NULL,
  fix = NULL,
  weights = NULL,
  family = "gaussian",
  model = "logistic",
  truncated = TRUE,
  link = NULL,
  b1type = c("asymptote", "rate"),
  method = "lm",
  gsl = NULL,
  gsu = NULL,
  gsn = NULL
)
```

## Arguments

- form:

  Formula, with the names of the response and the time variable, e.g.
  `weight ~ day`

- fform:

  Formula, list of formulas, or NULL. Fixed effects for the model
  parameters: a formula of the form b1 + b2 ~ Treatment, or a list of
  such formulas, e.g. list(b1 ~ Treatment, b2 ~ Treatment + Site).
  Parameters not listed have a single intercept (~ 1). A list of five
  right-hand side formulas, one for each of b1-b5 (X1 = ~ 1, ..., X5 = ~
  1), is also accepted.

- data:

  Data frame. Observations with missing values in the variables of the
  model or in the weights are left out.

- start:

  Numeric vector or NULL. Starting values for each model parameter. Only
  NA entries are replaced.

- fix:

  Logical vector or NULL. A vector with the same length as start. If
  TRUE, a model coefficient is fixed to the starting value. The
  threshold b4 and the asymmetry b5 are fixed by default, at 0 and 1
  unless given in start; fix=c(b4=FALSE) or fix=c(b5=FALSE) estimates
  them. Without truncation, b4 is not used and always fixed, at 0 unless
  given in start.

- weights:

  Numeric vector or NULL. Weights, one for each row of data.

- family:

  Character string. The distribution of the response, see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- model:

  Character string. The growth model: "logistic", "gompertz", "hill",
  "weibull", "normal", "gamma" or "linear", see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- truncated:

  Logical. If TRUE (default), the curve is truncated at a threshold b4,
  by default fixed at 0, see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- link:

  Named character vector or NULL. Links of the parameters b1-b5, see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- b1type:

  Character string. The meaning of b1, "asymptote" (default) or "rate"
  for the maximum growth rate, see
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- method:

  Character string. One of "lm" or "grid".

- gsl:

  Numeric vector or NULL. Lower bounds for the grid search, same length
  as start, on the scale of the coefficients (the rate 1/b3 for b3 with
  the inverse link); NA (or NULL) for the default range.

- gsu:

  Numeric vector or NULL. Upper bounds for the grid search, same length
  as start; NA (or NULL) for the default range.

- gsn:

  Integer, a single number or a vector, or NULL. Number of grid points
  for each parameter, a single number or a vector with the same length
  as start (default 10). A parameter with gsn = 1 stays at its "lm"
  starting value.

## Value

A named numeric vector of starting values, one for each model
coefficient.

## Details

Method "lm" linearises the growth curve: for a given asymptote b1, the
mean response at each time t, relative to b1, gives Finv(y/b1) = z(t),
with the inverse Finv of the function F of the model (qlogis for
logistic and hill, -log(-log(q)) for gompertz, log(-log(1 - q)) for
weibull) and the transformed time z(t) (see
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)),
which is linear in the time t (logistic, gompertz) or in log(t) (hill,
weibull). A linear regression of Finv(y/b1) on t or log(t) gives b2 and
b3 from its intercept and slope, e.g. b3 as the inverse of the slope and
b2 as the time at which the line crosses zero for the logistic model;
with b2 or b3 given in `start`, only the other one is estimated. The
weibull shape b3 starts at least at 1.1. The observed responses often do
not reach the asymptote yet, so b1 is profiled: for asymptotes from just
above the largest mean response to four times it, the linearised fit
gives b2 and b3, b1 is refitted by least squares for this curve, and the
asymptote with the smallest residual sum of squares is used. Mean
responses below 1\\ linearised curve and are left out. With covariates
in `fform`, the pooled starting value of each parameter is projected
onto its design matrix, giving the intercept or each group the pooled
value and contrasts a value of zero. The threshold b4, fixed at 0 by
default, starts when estimated where the curve without truncation
reaches 1\\ time with a positive count, and for the lognormal and beta
families before all times, as positive responses are impossible before
the threshold. The asymmetry b5 starts at 1, the symmetric curve, with
which b2 and b3 are linearised unless b5 is given in start; b4 without
truncation starts at 0.

Method "grid" evaluates the negative log-likelihood (without random
effects) at all combinations of a grid of parameter values and returns
the best grid point; parameters not included in the grid start at the
values of method "lm". Without bounds in `gsl` and `gsu`, default ranges
are used for b2 (the time range extended by a quarter at both ends,
positive with the log link) and b3 (1/4 to 4 times its "lm" starting
value, log-spaced; on the scale of its link). Default ranges are only
used for parameters with a single coefficient in `fform`; the asymptote
b1 is only included in the grid with user-specified bounds.

## Examples

``` r
gcTMBstart(height ~ DAS, data=FabaHeights25, model="gompertz")
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept) b5:(Intercept) 
#>      131.87454       27.92580       19.32536        0.00000        1.00000 
# with a given asymptote, from which b2 and b3 are linearised
gcTMBstart(height ~ DAS, data=FabaHeights25, start=c(b1=120))
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept) b5:(Intercept) 
#>      120.00000       33.47664       11.57527        0.00000        1.00000 
# best point of a grid of b2 and b3 values
gcTMBstart(height ~ DAS, data=FabaHeights25, method="grid", gsn=c(1, 10, 10, 1, 1))
#> b1:(Intercept) b2:(Intercept) b3:(Intercept) b4:(Intercept) b5:(Intercept) 
#>     118.820547      33.583333       9.643623       0.000000       1.000000 
```
