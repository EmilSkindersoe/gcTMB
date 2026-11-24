# Convergence of a gcTMB fit

Whether the fit converged to a maximum of the likelihood. A fit that
failed with an error and was returned by
`control=list(onerror="warning")` (see
[`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md))
has not converged, with its error as the reason. Otherwise: the negative
log-likelihood is finite, the optimiser reports convergence, the Hessian
of the parameters is positive definite (not checked for fits with
`control=list(sdreport=FALSE)`, see
[`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md)),
and the largest absolute gradient is below `gradtol`. By default, fits
with standard deviations of the random effects close to zero count as
converged; with `singular=TRUE` such singular fits do not (a standard
deviation below `sdtol`, or with a standard error of its log that is not
finite or above 100, as in `summary`). Note that the Hessian of a
singular fit is often not positive definite.

## Usage

``` r
converged(
  object,
  gradtol = 0.001,
  singular = FALSE,
  sdtol = 1e-08,
  lltol = 1e-06
)
```

## Arguments

- object:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- gradtol:

  Numeric. The largest acceptable absolute gradient of the negative
  log-likelihood, projected on the bounds of the coefficients.

- singular:

  Logical. If TRUE, singular fits, with a standard deviation of the
  random effects close to zero, count as not converged; FALSE (default)
  accepts them.

- sdtol:

  Numeric. With singular=TRUE, standard deviations of the random effects
  below sdtol make the fit singular.

- lltol:

  Numeric. For a fit restarted by `control=list(restart=TRUE)`: if the
  restart improved the log-likelihood by less than lltol, the optimiser
  message and the gradient are not checked.

## Value

TRUE or FALSE, with the attribute `reasons`, a character vector of the
failed checks (empty for a converged fit). NA for a model with
parameters set by
[`update.gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/update.gcTMB.md).

## Details

A fit with `control=list(restart=TRUE)` (see
[`gcTMBcontrol`](https://emilskindersoe.github.io/gcTMB/reference/gcTMBcontrol.md))
that was optimised a second time from its estimates, and whose
log-likelihood the restart improved by less than `lltol`, is at the
maximum: the convergence message of the optimiser and the gradient are
then not checked. With few groups and several correlated random effects,
the likelihood is often too flat in some directions for nlminb to report
convergence or reach a gradient below gradtol, although the maximum is
reached.

## Examples

``` r
fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
converged(fit)
#> [1] TRUE
```
