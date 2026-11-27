# Simulation

This vignette covers simulation studies: simulating data from a fitted
model, refitting it, and comparing the estimates with the true
parameters, including the autocorrelation parameter and the variance
function of the residuals. The correlation structures themselves,
temporal and spatial, are described in the vignette *Correlation
structures*.

## The sorghum and maize data

The examples use the biomass of `SorghumMaize` (see
[`?SorghumMaize`](https://emilskindersoe.github.io/gcTMB/reference/SorghumMaize.md)),
with the lognormal beta growth function of the vignette *Correlation
structures*: each combination of crop, input level and block is a plot,
and the zero biomass at sowing is left out.

``` r

library(gcTMB)
library(ggplot2)
SorghumMaize$plot <- interaction(SorghumMaize$Crop, SorghumMaize$Input, SorghumMaize$Block, drop = TRUE)
sm <- subset(SorghumMaize, Biomass > 0)
ggplot(SorghumMaize, aes(x = DOY, y = Biomass, group = plot, colour = Input)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Crop) +
  labs(x = "Day of the year", y = "Biomass (Mg/ha)", colour = "Input") +
  theme_minimal()
```

![](simulation_files/figure-html/unnamed-chunk-1-1.png)

## The model with autocorrelation

The residuals of a plot have the continuous-time AR(1) correlation
$`\rho^{dt}`$ between samplings $`dt`$ days apart,
`corCAR1(form = ~ DOY | plot)` (see the vignette *Correlation
structures*):

``` r

fitrec <- gcTMB(Biomass ~ DOY, fform = b1 + b2 ~ Crop * Input, rform = b1 ~ 1 | plot, data = sm,
                model = "beta", start = c(b4 = 141), family = "lognormal",
                correlation = corCAR1(value = 0.9, form = ~ DOY | plot))
corpars(fitrec)
#>        rho 
#> 0.09005289
```

The estimate of $`\rho`$ is a correlation of practically zero between
consecutive samplings. Below, simulations with a true autocorrelation
show that the model finds it when it is there.

## Simulation

A simulation study has four steps: fit a model to data, simulate new
data sets from the fitted model, refit the model to the simulated data,
and compare the estimates with the true parameters, the estimates of the
original fit.

### The true model

We fit the lognormal beta growth function to the biomass, with an
asymptote and a time of the fastest growth for each crop, and a random
asymptote of each plot. The start of growth is fixed at sowing,
$`b_4 = 141`$; the beta curve has no asymmetry $`b_5`$.

``` r

fitsim <- gcTMB(Biomass ~ DOY, fform = b1 + b2 ~ Crop, rform = b1 ~ 1 | plot, data = sm,
                model = "beta", start = c(b4 = 141), family = "lognormal")
fitsim
#> Call:
#> gcTMB(form = Biomass ~ DOY, fform = b1 + b2 ~ Crop, rform = b1 ~ 
#>     1 | plot, data = sm, start = c(b4 = 141), family = "lognormal", 
#>     model = "beta")
#> 
#> Model:  beta
#> Family: lognormal, identity link for all parameters
#> 
#> Coefficients:
#>       b1:(Intercept)          b1:CropMaize  b1:CropSweet sorghum  
#>              36.2386              -12.5753                0.6892  
#>       b2:(Intercept)          b2:CropMaize  b2:CropSweet sorghum  
#>             236.9786               -5.0948                2.5446  
#>       b3:(Intercept)  
#>             278.2859  
#> Fixed: b4:(Intercept) = 141
#> 
#> Random effects:
#>  Groups Name Std.Dev.
#>  plot   b1   3.58    
#> 
#> Residual standard deviation of log(y): 0.3006
converged(fitsim)
#> [1] TRUE
```

The fitted model is the true model of the simulation. Its parameters
are:

``` r

fixef(fitsim)
#>       b1:(Intercept)         b1:CropMaize b1:CropSweet sorghum 
#>           36.2385683          -12.5752974            0.6891833 
#>       b2:(Intercept)         b2:CropMaize b2:CropSweet sorghum 
#>          236.9786410           -5.0947681            2.5446297 
#>       b3:(Intercept)       b4:(Intercept) 
#>          278.2858609          141.0000000
VarCorr(fitsim)
#> $plot
#>         b1
#> b1 12.8291
#> attr(,"stddev")
#>       b1 
#> 3.581773 
#> attr(,"correlation")
#>    b1
#> b1  1
#> 
#> attr(,"sc")
#> [1] 0.3005549
sigma(fitsim)
#> [1] 0.3005549
```

Other true models can be made from it with
[`update()`](https://rdrr.io/r/stats/update.html),
e.g. `update(fitsim, fixef = c("b1:(Intercept)" = 30))` for a larger
asymptote of fiber sorghum, `update(fitsim, VarCorr = ...)` for another
variance of the random effects, or `update(fit, corpar = ...)` for
another autocorrelation, as below.

### Simulating new plots

We simulate 20 new plots of each crop, each sampled every 10 days from
day 150 to day 300. Each plot in `newdata` is a new curve, with its own
random effect drawn from the fitted variance.

``` r

simdata <- expand.grid(DOY = seq(150, 300, by = 10), plot = paste("plot", 1:60, sep = ""))
simdata$Crop <- rep(levels(SorghumMaize$Crop), each = 20)[as.integer(simdata$plot)]
simdata$Biomass <- simulate(fitsim, nsim = 1, seed = 1, newdata = simdata)[[1]]
head(simdata)
#>   DOY  plot          Crop    Biomass
#> 1 150 plot1 Fiber sorghum 0.03351688
#> 2 160 plot1 Fiber sorghum 0.18146333
#> 3 170 plot1 Fiber sorghum 0.86335123
#> 4 180 plot1 Fiber sorghum 1.76762994
#> 5 190 plot1 Fiber sorghum 2.77908489
#> 6 200 plot1 Fiber sorghum 6.29144647
```

### Refitting the model

The same model is fitted to the simulated data.

``` r

simfit <- gcTMB(Biomass ~ DOY, fform = b1 + b2 ~ Crop, rform = b1 ~ 1 | plot, data = simdata,
                model = "beta", start = c(b4 = 141), family = "lognormal")
converged(simfit)
#> [1] FALSE
#> attr(,"reasons")
#> [1] "the largest absolute gradient is 0.0021"
```

The 95% parametric bootstrap confidence intervals of the coefficients
and the variance parameters are compared with the true values, the
estimates of the original fit. `confint` reports the standard deviation
of the random effects and the residual standard deviation on their
natural scale, while `coef` has them on the log scale, so the true
values and the estimates are transformed in the same way.

``` r

ci <- confint(simfit, method = "boot", nsim = 50, seed = 1)
#> Computing bootstrap confidence intervals (50 simulations) ...
natural <- function(p){
  lg <- grepl("^log_", names(p))
  p[lg] <- exp(p[lg])
  p
}
truth <- natural(coef(fitsim))
tab <- data.frame(true = truth, estimate = natural(coef(simfit)), ci,
                  covered = truth >= ci[, 1] & truth <= ci[, 2],
                  check.names = FALSE)
# the standard deviations are on their natural scale: drop "log_" from their names
rownames(tab) <- sub("^log_", "", rownames(tab))
tab
#>                             true    estimate       2.5 %      97.5 % covered
#> b1:(Intercept)        36.2385683  36.6772959  34.4540215  39.3830376    TRUE
#> b1:CropMaize         -12.5752974 -13.7901245 -17.3928779 -10.6321687    TRUE
#> b1:CropSweet sorghum   0.6891833   0.6887689  -1.6154835   3.6212082    TRUE
#> b2:(Intercept)       236.9786410 236.8808074 235.7555142 238.3891990    TRUE
#> b2:CropMaize          -5.0947681  -5.0449775  -5.8784663  -4.0220119    TRUE
#> b2:CropSweet sorghum   2.5446297   2.8495673   2.0070914   3.5961009    TRUE
#> b3:(Intercept)       278.2858609 278.6154339 277.1445592 280.2344654    TRUE
#> sd:b1|plot             3.5817730   2.8538639   1.7905614   3.4507593   FALSE
#> sigma                  0.3005549   0.3154058   0.3011982   0.3311471   FALSE
```

The simulated biomass and the estimated curves of the 60 plots:

``` r

grid <- expand.grid(DOY = seq(141, 300, length.out = 200), plot = unique(simdata$plot))
grid$Crop <- rep(levels(SorghumMaize$Crop), each = 20)[as.integer(grid$plot)]
grid$Biomass <- predict(simfit, newdata = grid)
ggplot(simdata, aes(DOY, Biomass, group = plot)) +
  geom_line(data = grid, colour = "grey50", linewidth = 0.3) +
  geom_point(size = 0.8) +
  facet_wrap(~ Crop) +
  labs(x = "Day of the year", y = "Biomass (Mg/ha)") +
  theme_bw()
```

![](simulation_files/figure-html/sim-plot-1.png)

The simulated biomass is positive, and varies with the size of the
plants, as the observed biomass does: the lognormal residual standard
deviation, 0.3 on the log scale, is a coefficient of variation of about
30%. A simulation study reproduces the fitted model, not the process
that generated the data, so its conclusions only carry over to the real
data as far as the model describes them; with the gaussian family, about
a tenth of the simulated biomass values would be negative.

In a full simulation study, these steps are repeated for many simulated
data sets, e.g. to estimate the coverage of the confidence intervals or
to compare growth models with
[`criteria()`](https://emilskindersoe.github.io/gcTMB/reference/criteria.md),
as below for the autocorrelation.

### Recovering the autocorrelation

With a correlation structure,
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws the
residuals of each group with the autocorrelation of the fit: within a
plot in `newdata`, ordered in time, each residual given the previous
one, $`dt`$ earlier, is normal with mean $`\rho^{dt}`$ times the
previous residual and variance $`\sigma^2 (1 - \rho^{2 dt})`$.
[`update()`](https://rdrr.io/r/stats/update.html) is designed for
simulation studies like this one: it sets some of the parameters of a
fitted model to chosen values and keeps the others, without refitting,
so that the result can serve as the true model of the simulation. Its
arguments are `fixef` for the coefficients, `VarCorr` for the covariance
matrices of the random effects, `sigma` for the residual standard
deviation, `rho` for the autocorrelation, and `varpar` for the
parameters of a variance function.

Here the true model is the lognormal fit `fitrec` from above, with crop
and input effects, random asymptotes of the plots and continuous-time
AR(1) residuals of $`\log(y)`$. Its autocorrelation, $`\rho = 0.09`$ per
day, is set to $`\rho = 0.8`$ with
[`update()`](https://rdrr.io/r/stats/update.html), a correlation of
$`0.8^{10} \approx 0.11`$ between samplings 10 days apart, and all other
parameters keep their estimates:

``` r

rho0 <- 0.8
fitcar <- update(fitrec, corpar = rho0)
```

[`update()`](https://rdrr.io/r/stats/update.html) changes $`\rho`$ and
nothing else, so the other parameters are those estimated together with
$`\rho = 0.09`$. Alternatively,
[`refit()`](https://emilskindersoe.github.io/gcTMB/reference/refit.md)
fits the model again with $`\rho`$ fixed at 0.8,
`corCAR1(value = 0.8, fixed = TRUE)`, and estimates the other parameters
given $`\rho = 0.8`$. This is the way to go when the true model should
fit the data as well as possible with the chosen $`\rho`$, and it has
the higher log-likelihood;
[`update()`](https://rdrr.io/r/stats/update.html) is the way to go when
only $`\rho`$ should differ from the fitted model, e.g. to compare
simulations with different autocorrelations and otherwise the same
parameters.

``` r

fitfix <- refit(fitrec, correlation = corCAR1(value = rho0, form = ~ DOY | plot, fixed = TRUE))
pars <- function(f) c(fixef(f)[c("b1:(Intercept)", "b2:(Intercept)", "b3:(Intercept)")],
                      sigma = sigma(f), corpars(f), logLik = as.numeric(logLik(f)))
rbind(estimated = pars(fitrec), updated = pars(fitcar), refitted = pars(fitfix))
#>           b1:(Intercept) b2:(Intercept) b3:(Intercept)     sigma        rho
#> estimated       32.67956       237.4210       278.3070 0.2904354 0.09005289
#> updated         32.67956       237.4210       278.3070 0.2904354 0.80000000
#> refitted        32.60727       237.4475       278.3419 0.2938293 0.80000000
#>              logLik
#> estimated -466.2650
#> updated   -468.4116
#> refitted  -468.3101
```

We simulate from the updated model: 50 data sets of 10 new plots of each
combination of crop and input level, each sampled every 10 days from day
150 to day 300. We fit the model to each data set, now with $`\rho`$
estimated. [`confint()`](https://rdrr.io/r/stats/confint.html) computes
the Wald interval of $`\rho`$ on the logit scale and transforms it back.

``` r

design <- expand.grid(DOY = seq(150, 300, by = 10), rep = 1:10, Crop = levels(SorghumMaize$Crop),
                      Input = levels(SorghumMaize$Input))
design$plot <- interaction(design$rep, design$Crop, design$Input, drop = TRUE)
sims <- simulate(fitcar, nsim = 50, seed = 1, newdata = design)
carsim <- do.call(rbind, lapply(sims, function(y){
  dat <- cbind(design, Biomass = y)
  f1 <- gcTMB(Biomass ~ DOY, fform = b1 + b2 ~ Crop * Input, rform = b1 ~ 1 | plot, data = dat,
              model = "beta", start = c(b4 = 141), family = "lognormal",
              correlation = corCAR1(value = 0.9, form = ~ DOY | plot))
  ci <- confint(f1, parm = "cor_par:CAR1", method = "wald")
  data.frame(rho = corpars(f1)[["rho"]], lower = ci[1], upper = ci[2],
             covered = ci[1] <= rho0 & rho0 <= ci[2])
}))
```

The means of the estimates of $`\rho`$ and of the bounds of their 95%
Wald intervals, and the coverage of the intervals, against the true
value 0.8:

``` r

colMeans(carsim[, c("rho", "lower", "upper", "covered")])
#>      rho    lower    upper  covered 
#> 0.774071      NaN      NaN       NA
```

``` r

ggplot(carsim, aes(x = rho)) +
  geom_histogram(bins = 10, fill = "grey70", colour = "grey30") +
  geom_vline(xintercept = rho0, colour = "red") +
  labs(x = expression(paste("Estimate of ", rho)), y = "Simulations") +
  theme_minimal()
```

![](simulation_files/figure-html/sim-carplot-1.png)

### Recovering a variance function

A variance function (see
[`?varFunc`](https://emilskindersoe.github.io/gcTMB/reference/varFunc.md)
and the vignette *Variance functions*) makes the residual standard
deviation depend on a covariate, here on the fitted curve: with
`variance = varPower()`, the standard deviation of $`\log(y)`$ is
$`\sigma f(t)^\delta`$. A negative power $`\delta`$ makes the small
plants relatively more variable than the large ones. We fit it to the
biomass:

``` r

fitvar <- gcTMB(Biomass ~ DOY, fform = b1 + b2 ~ Crop, rform = b1 ~ 1 | plot, data = sm,
                model = "beta", start = c(b4 = 141), family = "lognormal", variance = varPower())
varpars(fitvar)
#>      power 
#> -0.2255527
AIC(fitsim, fitvar)
#>        df      AIC
#> fitsim  9 984.9414
#> fitvar 10 970.1678
```

The power is negative because it acts on the scale of $`\log(y)`$, where
the standard deviation is about the coefficient of variation: the
relative deviations of the biomass from the curve are larger for small
plants than for large ones. On the original scale, the standard
deviation still grows with the biomass, about as $`f(t)^{1 + \delta}`$;
the same model with the gaussian family estimates a power of about 0.8.

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws the
residuals with the standard deviation of the variance function at each
simulated curve. The true model sets the power to $`\delta = -0.3`$ with
[`update()`](https://rdrr.io/r/stats/update.html) and keeps the other
estimates:

``` r

delta0 <- -0.3
fitvartrue <- update(fitvar, varpar = c(power = delta0))
```

We simulate 50 data sets of 10 new plots of each crop, each sampled
every 10 days from day 150 to day 300, and fit the model to each, with
the power estimated:

``` r

design <- expand.grid(DOY = seq(150, 300, by = 10), rep = 1:10, Crop = levels(SorghumMaize$Crop))
design$plot <- interaction(design$rep, design$Crop, drop = TRUE)
sims <- simulate(fitvartrue, nsim = 50, seed = 1, newdata = design)
varsim <- do.call(rbind, lapply(sims, function(y){
  dat <- cbind(design, Biomass = y)
  f1 <- gcTMB(Biomass ~ DOY, fform = b1 + b2 ~ Crop, rform = b1 ~ 1 | plot, data = dat,
              model = "beta", start = c(b4 = 141), family = "lognormal", variance = varPower())
  ci <- confint(f1, parm = "var_par:power", method = "wald")
  data.frame(power = varpars(f1)[["power"]], lower = ci[1], upper = ci[2],
             covered = ci[1] <= delta0 & delta0 <= ci[2])
}))
```

The means of the estimates of $`\delta`$ and of the bounds of their 95%
Wald intervals, and the coverage of the intervals, against the true
value $`-0.3`$:

``` r

colMeans(varsim[, c("power", "lower", "upper", "covered")])
#>      power      lower      upper    covered 
#> -0.3004623 -0.3325021 -0.2684226  0.9200000
```

``` r

ggplot(varsim, aes(x = power)) +
  geom_histogram(bins = 10, fill = "grey70", colour = "grey30") +
  geom_vline(xintercept = delta0, colour = "red") +
  labs(x = expression(paste("Estimate of ", delta)), y = "Simulations") +
  theme_minimal()
```

![](simulation_files/figure-html/sim-varplot-1.png)
