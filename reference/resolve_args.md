# Concrete arguments for one design

Replaces each
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html) argument
by the design's value (looked up by id) and applies scenario overrides.

## Usage

``` r
resolve_args(object, values = list(), scenario = list())
```

## Arguments

- object:

  A simulation spec.

- values:

  Named list of tuned values, keyed by parameter id.

- scenario:

  Named list of scenario values (may be empty).

## Value

A named list of arguments.
