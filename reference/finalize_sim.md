# Finalize a spec with chosen values

Like `finalize_model()`: replaces every
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html) argument
with the value in `parameters` (a one-row tibble, for example from
[`select_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md)).

## Usage

``` r
finalize_sim(x, parameters)
```

## Arguments

- x:

  A simulation spec.

- parameters:

  A one-row data frame or named list of values keyed by id.

## Value

A spec with no tuned arguments.
