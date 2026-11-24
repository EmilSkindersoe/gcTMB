# Heights of faba bean plants in a field trial

Repeated height measurements of 320 faba bean plants in a field trial
(2026): five plants in each of 64 plots, with four genotypes on two soil
types, irrigated or not, measured on 14 days between 28 and 122 days
after sowing. The Loamy plots were not measured 77 days after sowing,
and five single measurements are missing, so the data have 4315 rows
instead of 4480.

## Usage

``` r
FabaHeights26
```

## Format

A data frame with 4315 rows and 7 variables:

- DAS:

  days after sowing

- Genotype:

  genotype: Lynx, Melodie, Navarra or Stella

- Plot_ID:

  identifier of the plot, 1-64

- Plant_ID:

  identifier of the plant, a factor unique to each plant: plot and plant
  number within the plot, e.g. "31_5"

- Soil:

  soil type: Loamy or Sandy

- Treatment:

  Irrigated or NotIrrigated

- Height:

  plant height in cm

The raw data and the script that converts them are in the directory
data-raw of the source package.

## Examples

``` r
fit <- gcTMB(Height ~ DAS, fform=b1 ~ Soil * Treatment, data=FabaHeights26)
summary(fit)
#> Call:
#> gcTMB(form = Height ~ DAS, fform = b1 ~ Soil * Treatment, data = FabaHeights26)
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 4315
#> Log-likelihood: -15267 (df = 7), AIC: 30548, BIC: 30593
#> 
#> Coefficients:
#>                                    Estimate Std. Error
#> b1:(Intercept)                       81.368      0.377
#> b1:SoilSandy                        -10.309      0.483
#> b1:TreatmentNotIrrigated             -6.871      0.496
#> b1:SoilSandy:TreatmentNotIrrigated  -12.455      0.683
#> b2:(Intercept)                       58.746      0.139
#> b3:(Intercept)                        8.119      0.123
#> Tests of the differences between curves, with p-values and confidence intervals: partest().
#> 
#> Residual standard deviation: 8.33 on 4309 degrees of freedom
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 6.3e-07
```
