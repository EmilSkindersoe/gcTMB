# Set control parameters for gcTMB

Set control parameters for gcTMB

## Usage

``` r
gcTMBcontrol(
  optimiser = "nlminb",
  iter.max = 1000,
  eval.max = 2000,
  smallsample = TRUE,
  sdreport = TRUE,
  restart = FALSE
)
```

## Arguments

- optimiser:

  Character string. The optimiser: only "nlminb"
  ([`nlminb`](https://rdrr.io/r/stats/nlminb.html), with the gradients
  of TMB) is available, for the fit, the profiles and the bootstrap.

- iter.max:

  Integer. Maximum number of iterations of the optimiser, see
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- eval.max:

  Integer. Maximum number of evaluations of the objective function, see
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- smallsample:

  Logical. If TRUE (default), gaussian and lognormal models without
  random effects and without autocorrelation use the unbiased residual
  variance with n - p residual degrees of freedom and t quantiles, as lm
  and nls; if FALSE, the maximum likelihood estimate and normal
  quantiles. See
  [`gcTMB`](https://emilskindersoe.github.io/gcTMB/reference/gcTMB.md).

- sdreport:

  Logical. If TRUE (default), the fit computes the covariance matrix of
  the estimates (TMB's sdreport), for standard errors, Wald intervals
  and the standard errors of predictions. FALSE skips it, which saves
  about a quarter of the time of a fit with random effects, e.g. in
  simulation studies that only need the estimates and the
  log-likelihood: the standard errors, `vcov` and the Wald intervals are
  then NA, and
  [`converged`](https://emilskindersoe.github.io/gcTMB/reference/converged.md)
  does not check the Hessian. Profile and bootstrap intervals are still
  available.

- restart:

  Logical. If TRUE, a fit whose optimiser does not report convergence,
  or whose largest absolute gradient is above 1e-3, is optimised once
  more from its estimates, and keeps the better of the two optima. With
  few groups and several correlated random effects, the likelihood is
  often so flat in some directions that nlminb reports "singular
  convergence" or "false convergence", or stops with a gradient slightly
  above the tolerance, at the maximum: the restart then hardly changes
  the log-likelihood. The improvement of the log-likelihood by the
  restart is kept in the fit (`fit$restart$improvement`), and
  [`converged`](https://emilskindersoe.github.io/gcTMB/reference/converged.md)
  accepts a fit that the restart improved by less than its argument
  lltol. FALSE (default) for no restart.

## Value

A list with the control parameters.

## Examples

``` r
gcTMBcontrol(iter.max=2000)
#> $optimiser
#> [1] "nlminb"
#> 
#> $iter.max
#> [1] 2000
#> 
#> $eval.max
#> [1] 2000
#> 
#> $smallsample
#> [1] TRUE
#> 
#> $sdreport
#> [1] TRUE
#> 
#> $restart
#> [1] FALSE
#> 
# small-sample inference switched off: maximum likelihood estimate of the residual sd
gcTMB(height ~ DAS, data=FabaHeights25, control=list(smallsample=FALSE))
#> Call:
#> gcTMB(form = height ~ DAS, data = FabaHeights25, control = list(smallsample = FALSE))
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>         120.79           30.73           12.88  
#> 
#> Residual standard deviation: 7.982
# estimates only, without standard errors
gcTMB(height ~ DAS, rform=b1 ~ 1 | plant_ID, data=FabaHeights25, control=list(sdreport=FALSE))
#> Call:
#> gcTMB(form = height ~ DAS, rform = b1 ~ 1 | plant_ID, data = FabaHeights25, 
#>     control = list(sdreport = FALSE))
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Coefficients:
#> b1:(Intercept)  b2:(Intercept)  b3:(Intercept)  
#>         120.56           30.61           12.89  
#> 
#> Random effects:
#>  Groups   Name Std.Dev.
#>  plant_ID b1   11.85   
#> 
#> Residual standard deviation: 3.127
# a second optimisation from the estimates when the first one does not converge
fit <- gcTMB(height ~ DAS, rform=b1 + b2 + b3 ~ 1 | plant_ID, data=FabaHeights25,
             control=list(restart=TRUE))
converged(fit)
#> [1] TRUE
```
