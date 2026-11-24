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

- New:

  number of new flowers

- Open:

  number of open flowers

- Old:

  number of old flowers

- Pods:

  number of pods

- N_flowers:

  number of flowers, the sum of New, Open and Old

The raw data and the script that converts them are in the directory
data-raw of the source package.

## Examples

``` r
fit <- gcTMB(Pods ~ DAS, fform=b1 ~ Genotype, data=FabaCounts23, family="poisson")
summary(fit)
#> Call:
#> gcTMB(form = Pods ~ DAS, fform = b1 ~ Genotype, data = FabaCounts23, 
#>     family = "poisson")
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: poisson, identity link for all parameters
#> 
#> Number of observations: 263
#> Log-likelihood: -445.7 (df = 8), AIC: 907.3, BIC: 935.9
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> b1:(Intercept)     18.17704    1.36539  13.313  < 2e-16 ***
#> b1:GenotypeCapri   -2.35787    1.60652  -1.468 0.142190    
#> b1:GenotypeILB938  -5.62794    1.54764  -3.636 0.000276 ***
#> b1:GenotypeLynx    -1.26468    1.62713  -0.777 0.437015    
#> b1:GenotypeMelodie  2.31500    1.82222   1.270 0.203931    
#> b1:GenotypeNavara  -2.48648    1.60413  -1.550 0.121129    
#> b2:(Intercept)     71.94586    0.22891 314.294  < 2e-16 ***
#> b3:(Intercept)      1.56828    0.09711  16.149  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Unadjusted Wald tests; tests of the differences between curves adjusted for their number: partest().
#> 
#> Optimiser: both X-convergence and relative convergence (5), maximum absolute gradient 1.1e-06
```
