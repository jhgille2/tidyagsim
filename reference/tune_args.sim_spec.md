# Tuned arguments of a spec

Tuned arguments of a spec

## Usage

``` r
# S3 method for class 'sim_spec'
tune_args(object, ...)
```

## Arguments

- object:

  A simulation spec.

- ...:

  Not used.

## Value

A tibble with one row per argument marked
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html): `name`,
`tunable`, `id`, `source` and `component`, as in tidymodels'
[`tune_args()`](https://generics.r-lib.org/reference/tune_args.html).
