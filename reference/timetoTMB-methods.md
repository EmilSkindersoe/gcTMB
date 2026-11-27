# timetoTMB Methods

timetoTMB Methods

## Usage

``` r
# S3 method for class 'timetoTMB'
print(x, n = NULL, ...)

# S3 method for class 'timetoTMB'
coef(object, ...)

# S3 method for class 'timetoTMB'
vcov(object, ...)
```

## Arguments

- x, object:

  A timetoTMB object, the result of
  [`timeto`](https://emilskindersoe.github.io/gcTMB/reference/timeto.md).

- n:

  Integer or NULL. The maximum number of times printed; NULL (default)
  for the option `gcTMB.print.max`, 100 unless set, e.g.
  `options(gcTMB.print.max = 200)`; `Inf` for all. When times are left
  out, a note at the end says how many, and how to show or extract all
  of them

- ...:

  Not used.

## Value

Estimated times to reach the levels of the asymptote, and their
covariance matrix
