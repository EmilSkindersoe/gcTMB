# Flower and pod counts of faba bean plants in a field trial

Repeated counts of flowers and pods of 96 faba bean plants in a field
trial (2026): three plants in each of 32 plots, with two genotypes on
two soil types, irrigated or not, counted on eight days between 59 and
84 days after sowing. The plants are plants 1-3 of the five of each plot
in
[`FabaHeights26`](https://emilskindersoe.github.io/gcTMB/reference/FabaHeights26.md),
with the same `Plant_ID`, so the two data sets can be merged. The counts
were recorded per node, each node only once: on the first date it had
flowers (or pods), as no more appear on a node afterwards. `Flowers` and
`Pods` are therefore the new flowers and pods since the previous date,
summed over the nodes of each plant, and `FlowersTotal` and `PodsTotal`
the running totals, which follow a growth curve. A plant without new
flowers or pods on a date, including a plant that was no longer
recorded, has 0 new flowers or pods. Flowering had begun before the
first count, so the flowers of the first date include those of earlier
days, and `FlowersTotal` does not start at 0.

## Usage

``` r
FabaFlowers26
```

## Format

A data frame with 768 rows and 10 variables:

- DAS:

  days after sowing

- Genotype:

  genotype: Navarra or Stella

- Plot_ID:

  identifier of the plot

- Plant_ID:

  identifier of the plant, a factor unique to each plant: plot and plant
  number within the plot, e.g. "31_3", as in `FabaHeights26`

- Soil:

  soil type: Loamy or Sandy

- Treatment:

  Irrigated or NotIrrigated

- Flowers:

  number of new flowers on the plant since the previous date

- Pods:

  number of new pods on the plant since the previous date

- FlowersTotal:

  number of flowers the plant has produced up to and including the date

- PodsTotal:

  number of pods the plant has set up to and including the date

The raw data and the script that converts them, `DIGIFABA_field_26.R`,
are in the directory data-raw of the source package.

## Examples

``` r
fit <- gcTMB(PodsTotal ~ DAS, fform=b1 + b2 ~ Soil * Treatment, rform=b1 ~ 1 | Plant_ID,
             data=FabaFlowers26)
summary(fit)
#> Call:
#> gcTMB(form = PodsTotal ~ DAS, fform = b1 + b2 ~ Soil * Treatment, 
#>     rform = b1 ~ 1 | Plant_ID, data = FabaFlowers26)
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 768, groups: Plant_ID 96
#> Log-likelihood: -1944 (df = 11), AIC: 3910, BIC: 3961
#> 
#> Coefficients:
#>                                    Estimate Std. Error
#> b1:(Intercept)                      27.1292      1.646
#> b1:SoilSandy                        -9.8595      2.303
#> b1:TreatmentNotIrrigated           -12.0468      2.302
#> b1:SoilSandy:TreatmentNotIrrigated   2.1591      3.245
#> b2:(Intercept)                      74.8665      0.165
#> b2:SoilSandy                        -0.6645      0.271
#> b2:TreatmentNotIrrigated            -1.2939      0.290
#> b2:SoilSandy:TreatmentNotIrrigated  -0.4712      0.605
#> b3:(Intercept)                       2.1007      0.085
#> Tests of the differences between curves, with p-values and confidence intervals: partest().
#> 
#> Random effects:
#>  Groups   Name Std.Dev.
#>  Plant_ID b1   7.754   
#> 
#> Residual standard deviation: 2.468
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 0.00032
timeto(fit, level=50, newdata=expand.grid(Soil=c("Loamy", "Sandy"),
                                          Treatment=c("Irrigated", "NotIrrigated")), random=FALSE)
#> Times to reach the percentages of the asymptote:
#> 
#>                                       Estimate Std.Error
#> 50%:Soil=Loamy,Treatment=Irrigated        74.9     0.165
#> 50%:Soil=Sandy,Treatment=Irrigated        74.2     0.230
#> 50%:Soil=Loamy,Treatment=NotIrrigated     73.6     0.252
#> 50%:Soil=Sandy,Treatment=NotIrrigated     72.4     0.485
```
