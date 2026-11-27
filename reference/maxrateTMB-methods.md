# maxrateTMB Methods

maxrateTMB Methods

## Usage

``` r
# S3 method for class 'maxrateTMB'
print(x, n = NULL, ...)

# S3 method for class 'maxrateTMB'
coef(object, ...)

# S3 method for class 'maxrateTMB'
vcov(object, ...)
```

## Arguments

- x, object:

  A maxrateTMB object, the result of
  [`maxrate`](https://emilskindersoe.github.io/gcTMB/reference/maxrate.md).

- n:

  Integer or NULL. The maximum number of estimates printed; NULL
  (default) for the option `gcTMB.print.max`, 100 unless set, e.g.
  `options(gcTMB.print.max = 200)`; `Inf` for all. When estimates are
  left out, a note at the end says how many, and how to show or extract
  all of them

- ...:

  Not used.

## Value

The estimates of the timing and value of the maximum growth rate, and
their covariance matrix
