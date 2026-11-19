# Likelihood profiles of a gcTMB model

Profiles the (marginal) log-likelihood for single parameters by fixing
the parameter at a sequence of values and re-optimising all other
parameters, using
[`tmbprofile`](https://rdrr.io/pkg/TMB/man/tmbprofile.html). With random
effects, the profiled likelihood is the Laplace approximation of the
marginal likelihood.

## Usage

``` r
# S3 method for class 'gcTMB'
profile(
  fitted,
  parm = NULL,
  level = 0.99,
  trace = FALSE,
  progress = interactive(),
  ...
)

# S3 method for class 'profile.gcTMB'
confint(object, parm = NULL, level = 0.95, ...)

# S3 method for class 'profile.gcTMB'
plot(x, level = 0.95, ...)
```

## Arguments

- fitted:

  A gcTMB object, the result of
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- parm:

  Character or integer vector, or NULL. Parameters to profile, names or
  positions in `coef(fitted)`. Default are the estimated coefficients
  b1-b5, without the variance parameters.

- level:

  Numeric between 0 and 1. The highest confidence level the profiles are
  computed for: each profile extends until the deviance from the maximum
  exceeds the cutoff (see Details). For confint and plot, the confidence
  level of the intervals and of the cutoff line.

- trace:

  Logical. Print the progress of
  [`tmbprofile`](https://rdrr.io/pkg/TMB/man/tmbprofile.html)

- progress:

  Logical. Show a progress bar over the profiled parameters

- ...:

  Further arguments passed to
  [`tmbprofile`](https://rdrr.io/pkg/TMB/man/tmbprofile.html), e.g.
  `ystep`; ignored by the confint and plot methods

- object, x:

  A profile.gcTMB object, the result of `profile`.

## Value

profile returns a profile.gcTMB object, a list with a data frame for
each parameter with the parameter values (`par`) and the deviance from
the maximum (`deviance`), i.e. twice the difference of the negative
log-likelihoods. confint returns a matrix of confidence intervals. plot
returns a ggplot object (package ggplot2) with a panel for each
parameter, showing the deviance, the cutoff for `level` (dashed), and
the estimate (dotted).

## Details

The profiles are computed on the scale of `coef`: log standard
deviations of the random effects (log_sd), parameters of their
correlation matrix (corr_par), and the log residual standard deviation
(log_sigma), log precision (log_phi), or log size (log_theta) or
dispersion (log_phi) of the negative binomial families. The confint
method transforms the intervals of the log standard deviations and log
precision to the natural scale; correlation parameters are only
transformed for terms with two random effects, as their parameter
determines the correlation only then. If the profile of a log standard
deviation does not reach the cutoff at low values, the likelihood is
flat towards zero and the lower limit of the standard deviation is 0.
Bounds given by `lower` and `upper` in
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md) are
used for the profiled coefficient, but not when re-optimising the other
coefficients.

The cutoff of the deviance is the quantile of the chi-square
distribution with one degree of freedom. For gaussian and lognormal
models without random effects and small-sample inference (see
[`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md)),
the cutoff of the coefficients b1-b5 is F-based instead, sum(w) log(1 +
F/(n - p)) with the quantile F of the F distribution with 1 and n - p
degrees of freedom, which gives the intervals of a t-test of each
coefficient, as profile intervals of nls.

## Examples

``` r
fit <- gcTMB(height ~ DAS, data=FabaHeights25)
pr <- profile(fit, progress=FALSE)
confint(pr)
#>                    2.5 %    97.5 %
#> b1:(Intercept) 114.54544 129.48091
#> b2:(Intercept)  29.43847  32.18230
#> b3:(Intercept)  11.21241  15.18937
plot(pr)
```
