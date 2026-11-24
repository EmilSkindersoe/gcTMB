# Flower and pod counts of faba bean plants

Repeated counts of flowers and pods of 22 faba bean plants of six
genotypes, all with the treatment "ww", on twelve days between 61 and 79
days after sowing (2023).

## Usage

``` r
FabaCounts23
```

## Format

A data frame with 263 rows and 9 variables:

- Tube_ID:

  identifier of the tube of the plant

- Genotype:

  genotype: Aurora, Capri, ILB938, Lynx, Melodie or Navara

- Treatment:

  treatment, "ww" for all plants

- DAS:

  days after sowing

- new:

  number of new flowers

- open:

  number of open flowers

- old:

  number of old flowers

- pods:

  number of pods

- n_flowers:

  number of flowers, the sum of new, open and old

The raw data and the script that converts them are in the directory
data-raw of the source package.

## Examples

``` r
fit <- gcTMB(pods ~ DAS, fform=b1 ~ Genotype, data=FabaCounts23, family="poisson")
summary(fit)
#> Call:
#> gcTMB(form = pods ~ DAS, fform = b1 ~ Genotype, data = FabaCounts23, 
#>     family = "poisson")
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: poisson, identity link for all parameters
#> 
#> Number of observations: 263
#> Log-likelihood: -445.7 (df = 8), AIC: 907.3, BIC: 935.9
#> 
#> Coefficients:
#>                    Estimate Std. Error
#> b1:(Intercept)       18.177      1.365
#> b1:GenotypeCapri     -2.358      1.607
#> b1:GenotypeILB938    -5.628      1.548
#> b1:GenotypeLynx      -1.265      1.627
#> b1:GenotypeMelodie    2.315      1.822
#> b1:GenotypeNavara    -2.486      1.604
#> b2:(Intercept)       71.946      0.229
#> b3:(Intercept)        1.568      0.097
#> Tests of the differences between curves, with p-values and confidence intervals: partest().
#> 
#> Optimiser: both X-convergence and relative convergence (5), maximum absolute gradient 1.1e-06
```
