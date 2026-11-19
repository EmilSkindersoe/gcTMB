# Timing and maximum growth rates

This vignette covers the times at which growth curves reach a level of
their asymptote,
[`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md),
and the maximum growth rates and their timing,
[`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md),
with treatment effects for the heights of faba bean plants and genotype
effects for counts of pods. The table of the inflection points and
maximum growth rates of the growth models is in the vignette *Growth
models and their parameters*.

The models of the vignettes *Growth models and their parameters* and
*Fixed and random effects*: a single logistic curve, a curve for each
treatment with random effects of the plants, and a Poisson model of the
pods with a curve for each genotype and a random asymptote of each
plant:

``` r

library(gcTMB)
library(ggplot2)
FabaHeights25$trt <- factor(FabaHeights25$Treatment)
trts <- data.frame(trt = levels(FabaHeights25$trt))
fit <- gcTMB(height ~ DAS, data = FabaHeights25, model = "logistic")
fitre <- gcTMB(height ~ DAS, fform = b1 + b2 ~ trt, rform = b1 + b2 ~ 1 | plant_ID, data = FabaHeights25)
fitp <- gcTMB(pods ~ DAS, fform = b1 + b2 ~ Genotype, rform = b1 ~ 1 | Tube_ID, data = FabaCounts23,
              family = "poisson", fix = c(b4 = FALSE))
genotypes <- data.frame(Genotype = sort(unique(FabaCounts23$Genotype)))
```

## Time to reach a level of the asymptote

[`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
gives the times at which a curve reaches a percentage of its asymptote,
$`b_2 + b_3 F^{-1}(p/100)`$, with standard errors by automatic
differentiation. With `type = "absolute"`, the times at which the curve
reaches given responses, here a height of 50cm.

``` r

tt <- timeto(fit, level = c(10, 50, 90))
tt
#> Times to reach the percentages of the asymptote:
#> 
#>     Estimate Std.Error
#> 10%     10.8     0.423
#> 50%     32.9     0.911
#> 90%     60.3     2.743
timeto(fit, level = 50, type = "absolute")
#> Times to reach the responses:
#> 
#>      Estimate Std.Error
#> y=50     28.8     0.336
```

[`confint()`](https://rdrr.io/r/stats/confint.html) gives profile
likelihood intervals by default, `method = "wald"` Wald intervals and
`method = "boot"` parametric bootstrap intervals.

``` r

confint(tt)
#>        2.5 %   97.5 %
#> 10% 10.02726 11.69745
#> 50% 31.32441 35.05932
#> 90% 55.61967 66.81695
confint(tt, method = "wald")
#>         2.5 %   97.5 %
#> 10%  9.997343 11.66661
#> 50% 31.109239 34.70278
#> 90% 54.875711 65.69606
```

## Maximum growth rate

[`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
gives the time of the maximum growth rate (`tmax`), the inflection point
of the curve, the maximum growth rate there (`rate`), and the proportion
of total growth at that time, $`f(t_m)/b_1`$ (`fraction`). The standard
errors are by the delta method.
[`confint()`](https://rdrr.io/r/stats/confint.html) gives profile
likelihood intervals for `tmax` and `rate` by default, and Wald
intervals for `fraction`; `method = "wald"` and `method = "boot"` as for
[`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).
For the logistic curve truncated at $`b_4 = 0`$, `tmax` is $`b_2`$;
`rate` is slightly above and `fraction` slightly below the values
$`b_1/(4 b_3)`$ and $`1/2`$ of the curve without truncation in the table
of the growth models, as the truncation shifts the curve down to zero at
$`t = 0`$ and rescales it.

``` r

mr <- maxrate(fit)
mr
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>          Estimate Std.Error
#> tmax       30.730   0.68220
#> rate        2.560   0.09415
#> fraction    0.454   0.00778
confint(mr)
#>               2.5 %     97.5 %
#> tmax     29.4384590 32.1825247
#> rate      2.3803996  2.7522120
#> fraction  0.4386425  0.4693352
```

## Timing and maximum growth rates of the treatments

[`timeto()`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md)
and
[`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
with `newdata` give the times and rates of the curves of each row, here
the population curves of the treatments (`random = FALSE`) of the model
with random effects of the plants:

``` r

timeto(fitre, level = 90, newdata = trts, random = FALSE)
#> Times to reach the percentages of the asymptote:
#> 
#>             Estimate Std.Error
#> 90%:trt=0.1     60.3      1.24
#> 90%:trt=0.3     58.9      1.21
#> 90%:trt=0.5     60.4      1.24
#> 90%:trt=0.6     60.9      1.25
#> 90%:trt=0.8     61.1      1.25
mrre <- maxrate(fitre, newdata = trts)
mrre
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>                  Estimate Std.Error
#> tmax:trt=0.1       30.817   1.14107
#> rate:trt=0.1        2.457   0.09484
#> fraction:trt=0.1    0.455   0.00446
#> tmax:trt=0.3       29.322   1.13761
#> rate:trt=0.3        2.441   0.09632
#> fraction:trt=0.3    0.449   0.00501
#> tmax:trt=0.5       30.995   1.13233
#> rate:trt=0.5        2.661   0.09567
#> fraction:trt=0.5    0.455   0.00436
#> tmax:trt=0.6       31.492   1.13281
#> rate:trt=0.6        2.671   0.09521
#> fraction:trt=0.6    0.457   0.00418
#> tmax:trt=0.8       31.720   1.13525
#> rate:trt=0.8        2.618   0.09471
#> fraction:trt=0.8    0.458   0.00412
```

### Treatment effects on the maximum growth rate

With `b1type = "rate"`, $`b_1`$ is the maximum growth rate of the curve
instead of its asymptote (see *The maximum growth rate as a parameter*
in the vignette *Growth models and their parameters*), so that
covariates on $`b_1`$ are effects on the maximum growth rate. The model
with the treatment effects on $`b_1`$ and $`b_2`$ has the same number of
parameters as `fitre`, but the random effects of the plants are now on
their maximum growth rates rather than their asymptotes, a slightly
different model. The log link gives the treatment effects as ratios of
the rates.

``` r

fitrr <- gcTMB(height ~ DAS, fform = b1 + b2 ~ trt, rform = b1 + b2 ~ 1 | plant_ID, data = FabaHeights25,
               b1type = "rate", link = c(b1 = "log"))
AIC(fitre, fitrr)
#>       df      AIC
#> fitre 15 1046.344
#> fitrr 15 1047.093
```

Whether the treatments have the same maximum growth rate is a likelihood
ratio test of the models with and without the treatment effects on
$`b_1`$, keeping those on $`b_2`$:

``` r

fitrr0 <- gcTMB(height ~ DAS, fform = b2 ~ trt, rform = b1 + b2 ~ 1 | plant_ID, data = FabaHeights25,
                b1type = "rate", link = c(b1 = "log"))
anova(fitrr0, fitrr)
#> Likelihood ratio tests of gcTMB models
#> 
#>        npar    AIC    BIC  logLik deviance    LR Df Pr(>LR)
#> fitrr0   11 1044.1 1080.4 -511.06   1022.1                 
#> fitrr    15 1047.1 1096.6 -508.55   1017.1 5.024  4  0.2848
```

[`pairwise()`](https://emilskindersoe.github.io/gcTMB/reference/pairwise.md)
compares the maximum growth rates of all pairs of treatments (see the
pods below for the details); with the log link, as ratios:

``` r

pairwise(fitrr, newdata = trts)
#> Pairwise differences of b1 (maximum growth rate) on the scale of its log link:
#> 
#>                   Estimate Std.Error z value p.value  lower  upper ratio
#> trt=0.1 - trt=0.3  0.00528    0.0523   0.101   1.000 -0.137 0.1480 1.005
#> trt=0.1 - trt=0.5 -0.07783    0.0523  -1.488   0.570 -0.220 0.0648 0.925
#> trt=0.1 - trt=0.6 -0.08482    0.0523  -1.621   0.485 -0.228 0.0579 0.919
#> trt=0.1 - trt=0.8 -0.06376    0.0523  -1.218   0.741 -0.207 0.0790 0.938
#> trt=0.3 - trt=0.5 -0.08310    0.0522  -1.591   0.504 -0.226 0.0593 0.920
#> trt=0.3 - trt=0.6 -0.09010    0.0522  -1.724   0.419 -0.233 0.0524 0.914
#> trt=0.3 - trt=0.8 -0.06904    0.0523  -1.321   0.679 -0.212 0.0735 0.933
#> trt=0.5 - trt=0.6 -0.00699    0.0522  -0.134   1.000 -0.149 0.1354 0.993
#> trt=0.5 - trt=0.8  0.01406    0.0522   0.269   0.999 -0.128 0.1565 1.014
#> trt=0.6 - trt=0.8  0.02106    0.0523   0.403   0.995 -0.121 0.1636 1.021
#>                   ratio.lower ratio.upper
#> trt=0.1 - trt=0.3       0.872        1.16
#> trt=0.1 - trt=0.5       0.802        1.07
#> trt=0.1 - trt=0.6       0.797        1.06
#> trt=0.1 - trt=0.8       0.813        1.08
#> trt=0.3 - trt=0.5       0.798        1.06
#> trt=0.3 - trt=0.6       0.792        1.05
#> trt=0.3 - trt=0.8       0.809        1.08
#> trt=0.5 - trt=0.6       0.861        1.15
#> trt=0.5 - trt=0.8       0.879        1.17
#> trt=0.6 - trt=0.8       0.886        1.18
#> 
#> Wald z tests; p-values and confidence limits adjusted by the single-step (Tukey) method, 95% confidence limits, and ratios.
```

The treatments do not differ clearly in their maximum growth rates.

## Timing and maximum growth rates of the genotypes

The same functions apply to count data. For the pods, the times at which
the population curve of each genotype reaches half of its asymptote, and
the maximum growth rates, the most pods gained per day, of the Poisson
model `fitp` with genotype effects on the asymptote:

``` r

timeto(fitp, level = 50, newdata = genotypes, random = FALSE)
#> Times to reach the percentages of the asymptote:
#> 
#>                      Estimate Std.Error
#> 50%:Genotype=Aurora      70.1     0.459
#> 50%:Genotype=Capri       69.3     0.440
#> 50%:Genotype=ILB938      69.2     0.519
#> 50%:Genotype=Lynx        72.8     0.393
#> 50%:Genotype=Melodie     73.3     0.390
#> 50%:Genotype=Navara      72.7     0.401
maxrate(fitp, newdata = genotypes)
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>                           Estimate Std.Error
#> tmax:Genotype=Aurora        70.053  0.468918
#> rate:Genotype=Aurora         2.789  0.488091
#> fraction:Genotype=Aurora     0.489  0.005257
#> tmax:Genotype=Capri         69.203  0.464850
#> rate:Genotype=Capri          2.181  0.400919
#> fraction:Genotype=Capri      0.478  0.009656
#> tmax:Genotype=ILB938        69.045  0.553914
#> rate:Genotype=ILB938         1.786  0.382493
#> fraction:Genotype=ILB938     0.476  0.011629
#> tmax:Genotype=Lynx          72.797  0.391550
#> rate:Genotype=Lynx           3.667  0.506615
#> fraction:Genotype=Lynx       0.499  0.000776
#> tmax:Genotype=Melodie       73.298  0.388664
#> rate:Genotype=Melodie        4.733  0.639974
#> fraction:Genotype=Melodie    0.499  0.000580
#> tmax:Genotype=Navara        72.654  0.399801
#> rate:Genotype=Navara         3.331  0.485199
#> fraction:Genotype=Navara     0.498  0.000860
```

### Genotype effects on the maximum growth rate

The genotype effects of `fitp` are on the asymptote. With
`b1type = "rate"`, $`b_1`$ is the maximum growth rate, so the same
formulas give genotype effects and random effects of the plants on the
maximum growth rate, the most pods gained per day; the log link
expresses the genotype effects as ratios of the rates. The random
effects are now on the rate of each plant rather than its asymptote,
which is a different model, here with a somewhat worse fit by AIC.

``` r

fitpr <- gcTMB(pods ~ DAS, fform = b1 + b2 ~ Genotype, rform = b1 ~ 1 | Tube_ID, data = FabaCounts23,
               family = "poisson", fix = c(b4 = FALSE), b1type = "rate", link = c(b1 = "log"))
summary(fitpr)
#> Call:
#> gcTMB(form = pods ~ DAS, fform = b1 + b2 ~ Genotype, rform = b1 ~ 
#>     1 | Tube_ID, data = FabaCounts23, fix = c(b4 = FALSE), family = "poisson", 
#>     link = c(b1 = "log"), b1type = "rate")
#> 
#> Model:  logistic, truncated at the threshold b4; b1 is the maximum growth rate
#> Family: poisson, log link for b1, identity link for the other parameters
#> 
#> Number of observations: 263, groups: Tube_ID 22
#> Log-likelihood: -356.4 (df = 15), AIC: 742.8, BIC: 796.3
#> 
#> Coefficients:
#>                    Estimate Std. Error
#> b1:(Intercept)       1.0042      0.191
#> b1:GenotypeCapri    -0.3002      0.232
#> b1:GenotypeILB938   -0.4572      0.231
#> b1:GenotypeLynx      0.2753      0.231
#> b1:GenotypeMelodie   0.5368      0.245
#> b1:GenotypeNavara    0.1767      0.232
#> b2:(Intercept)      70.0041      0.469
#> b2:GenotypeCapri    -0.6698      0.635
#> b2:GenotypeILB938   -1.0464      0.711
#> b2:GenotypeLynx      2.7625      0.553
#> b2:GenotypeMelodie   3.3139      0.556
#> b2:GenotypeNavara    2.6192      0.561
#> b3:(Intercept)       1.3350      0.148
#> b4:(Intercept)      65.0000      0.024
#> The coefficients of b1 are on the log scale.
#> 
#> Random effects:
#>  Groups  Name Std.Dev.
#>  Tube_ID b1   0.2643  
#> 
#> Pearson dispersion: 0.4938 (chi-square 122.5 on 248 df)
#> 
#> Optimiser: converged, with the threshold b4 at an observation time (a kink of the likelihood), maximum absolute gradient 0.28
AIC(fitp, fitpr)
#>       df      AIC
#> fitp  15 737.1423
#> fitpr 15 742.7615
```

Whether the genotypes have the same maximum growth rate is a likelihood
ratio test of the models with and without the genotype effects on
$`b_1`$, keeping those on $`b_2`$:

``` r

fitpr0 <- gcTMB(pods ~ DAS, fform = b2 ~ Genotype, rform = b1 ~ 1 | Tube_ID, data = FabaCounts23,
                family = "poisson", fix = c(b4 = FALSE), b1type = "rate", link = c(b1 = "log"))
anova(fitpr0, fitpr)
#> Likelihood ratio tests of gcTMB models
#> 
#>        npar    AIC    BIC  logLik deviance     LR Df Pr(>LR)  
#> fitpr0   10 742.32 778.04 -361.16   722.32                    
#> fitpr    15 742.76 796.34 -356.38   712.76 9.5615  5 0.08866 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

[`pairwise()`](https://emilskindersoe.github.io/gcTMB/reference/pairwise.md)
compares the maximum growth rates of all pairs of genotypes by Wald
tests. By default, the p-values and the confidence intervals are
adjusted for the number of pairs by the single-step method, the
generalisation of Tukey’s honestly significant differences to the
correlated estimates of a model, so that the intervals cover all the
differences simultaneously; `adjust = "holm"` and the other methods of
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) adjust the
p-values only. With the log link, the differences are log ratios of the
rates, and the ratios are given as well. The same works for the other
parameters, e.g. `parm = "b2"` for the time of the maximum growth rate,
which is $`b_2`$ for the logistic curve.

``` r

pairwise(fitpr, newdata = genotypes)
#> Pairwise differences of b1 (maximum growth rate) on the scale of its log link:
#> 
#>                                    Estimate Std.Error z value p.value  lower
#> Genotype=Aurora - Genotype=Capri     0.3002     0.232   1.296 0.78771 -0.360
#> Genotype=Aurora - Genotype=ILB938    0.4572     0.231   1.976 0.35585 -0.202
#> Genotype=Aurora - Genotype=Lynx     -0.2753     0.231  -1.192 0.84196 -0.933
#> Genotype=Aurora - Genotype=Melodie  -0.5368     0.245  -2.191 0.24139 -1.235
#> Genotype=Aurora - Genotype=Navara   -0.1767     0.232  -0.762 0.97484 -0.838
#> Genotype=Capri - Genotype=ILB938     0.1570     0.218   0.722 0.98021 -0.463
#> Genotype=Capri - Genotype=Lynx      -0.5754     0.217  -2.649 0.08514 -1.194
#> Genotype=Capri - Genotype=Melodie   -0.8370     0.232  -3.610 0.00443 -1.498
#> Genotype=Capri - Genotype=Navara    -0.4769     0.218  -2.186 0.24385 -1.099
#> Genotype=ILB938 - Genotype=Lynx     -0.7324     0.217  -3.376 0.00970 -1.351
#> Genotype=ILB938 - Genotype=Melodie  -0.9940     0.232  -4.288 0.00032 -1.655
#> Genotype=ILB938 - Genotype=Navara   -0.6339     0.218  -2.907 0.04186 -1.255
#> Genotype=Lynx - Genotype=Melodie    -0.2616     0.231  -1.130 0.87004 -0.921
#> Genotype=Lynx - Genotype=Navara      0.0985     0.218   0.453 0.99757 -0.521
#> Genotype=Melodie - Genotype=Navara   0.3601     0.232   1.549 0.63192 -0.302
#>                                      upper ratio ratio.lower ratio.upper
#> Genotype=Aurora - Genotype=Capri    0.9600 1.350       0.698       2.612
#> Genotype=Aurora - Genotype=ILB938   1.1165 1.580       0.817       3.054
#> Genotype=Aurora - Genotype=Lynx     0.3828 0.759       0.393       1.466
#> Genotype=Aurora - Genotype=Melodie  0.1613 0.585       0.291       1.175
#> Genotype=Aurora - Genotype=Navara   0.4842 0.838       0.433       1.623
#> Genotype=Capri - Genotype=ILB938    0.7769 1.170       0.629       2.175
#> Genotype=Capri - Genotype=Lynx      0.0434 0.562       0.303       1.044
#> Genotype=Capri - Genotype=Melodie  -0.1763 0.433       0.224       0.838
#> Genotype=Capri - Genotype=Navara    0.1448 0.621       0.333       1.156
#> Genotype=ILB938 - Genotype=Lynx    -0.1142 0.481       0.259       0.892
#> Genotype=ILB938 - Genotype=Melodie -0.3335 0.370       0.191       0.716
#> Genotype=ILB938 - Genotype=Navara  -0.0127 0.531       0.285       0.987
#> Genotype=Lynx - Genotype=Melodie    0.3978 0.770       0.398       1.489
#> Genotype=Lynx - Genotype=Navara     0.7185 1.104       0.594       2.051
#> Genotype=Melodie - Genotype=Navara  1.0223 1.433       0.739       2.780
#> 
#> Wald z tests; p-values and confidence limits adjusted by the single-step (Tukey) method, 95% confidence limits, and ratios.
```

The overall likelihood ratio test is borderline, while the pairwise
comparisons single out the largest differences: Melodie gains pods
faster than ILB938 and Capri, and ILB938 more slowly than Lynx and
Navara. The Wald tests rely on the normal approximation of the
estimates, which is rougher than the likelihood ratio test with only 22
plants.

[`predict()`](https://rdrr.io/r/stats/predict.html) with
`type = "parameters"` gives the asymptote of the population curves as
well, and
[`maxrate()`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md)
the rates, which are now $`b_1`$:

``` r

cbind(genotypes, predict(fitpr, newdata = genotypes, type = "parameters", random = FALSE))
#>   Genotype       b1       b2       b3 b4 asymptote
#> 1   Aurora 2.729828 70.00409 1.335018 65 14.241985
#> 2    Capri 2.021915 69.33432 1.335018 65 10.392846
#> 3   ILB938 1.728166 68.95769 1.335018 65  8.775831
#> 4     Lynx 3.594824 72.76662 1.335018 65 19.139689
#> 5  Melodie 4.669520 73.31796 1.335018 65 24.886596
#> 6   Navara 3.257544 72.62327 1.335018 65 17.338102
maxrate(fitpr, newdata = genotypes)
#> Time of the maximum growth rate (tmax), maximum growth rate (rate) and proportion of total growth at inflection (fraction):
#> 
#>                           Estimate Std.Error
#> tmax:Genotype=Aurora        70.004  0.469410
#> rate:Genotype=Aurora         2.730  0.521222
#> fraction:Genotype=Aurora     0.488  0.005415
#> tmax:Genotype=Capri         69.334  0.461599
#> rate:Genotype=Capri          2.022  0.350820
#> fraction:Genotype=Capri      0.481  0.008514
#> tmax:Genotype=ILB938        68.958  0.559069
#> rate:Genotype=ILB938         1.728  0.297596
#> fraction:Genotype=ILB938     0.474  0.012515
#> tmax:Genotype=Lynx          72.767  0.387529
#> rate:Genotype=Lynx           3.595  0.616811
#> fraction:Genotype=Lynx       0.499  0.000783
#> tmax:Genotype=Melodie       73.318  0.392049
#> rate:Genotype=Melodie        4.670  0.883139
#> fraction:Genotype=Melodie    0.499  0.000559
#> tmax:Genotype=Navara        72.623  0.396267
#> rate:Genotype=Navara         3.258  0.564670
#> fraction:Genotype=Navara     0.498  0.000868
```
