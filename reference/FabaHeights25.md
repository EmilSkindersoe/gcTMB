# Heights of faba bean plants

Repeated height measurements of 25 faba bean plants, five at each of
five treatment levels, measured on eight days between 14 and 61 days
after sowing (2025).

## Usage

``` r
FabaHeights25
```

## Format

A data frame with 200 rows and 7 variables:

- Tube_ID:

  identifier of the tube the plant grows in

- Treatment:

  treatment level: 0.1, 0.3, 0.5, 0.6 or 0.8

- Date:

  date of the measurement

- DAS:

  days after sowing

- height:

  plant height

- DASfact:

  days after sowing, the same values as DAS

- plant_ID:

  identifier of the plant, the same values as Tube_ID

The raw data and the script that converts them are in the directory
data-raw of the source package.

## Examples

``` r
fit <- gcTMB(height ~ DAS, rform=b1 + b2 ~ 1 | plant_ID, data=FabaHeights25)
summary(fit)
#> Call:
#> gcTMB(form = height ~ DAS, rform = b1 + b2 ~ 1 | plant_ID, data = FabaHeights25)
#> 
#> Model:  logistic, truncated at b4 = 0
#> Family: gaussian, identity link for all parameters
#> 
#> Number of observations: 200, groups: plant_ID 25
#> Log-likelihood: -513.7 (df = 7), AIC: 1041, BIC: 1064
#> 
#> Coefficients:
#>                Estimate Std. Error
#> b1:(Intercept)   121.09      2.257
#> b2:(Intercept)    30.87      0.542
#> b3:(Intercept)    12.84      0.250
#> 
#> Random effects:
#>  Groups   Name Std.Dev. Corr 
#>  plant_ID b1   10.237        
#>           b2    2.561   -0.12
#> 
#> Residual standard deviation: 2.074
#> 
#> Optimiser: relative convergence (4), maximum absolute gradient 8.4e-05
```
