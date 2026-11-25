# Paired comparison of two designs

The difference `b - a` in a metric on shared seeds (and scenarios), with
its standard error and the one-sided p-value that `b` is worse than `a`.
Pairing on seeds is what makes small differences between designs
detectable.

## Usage

``` r
compare_designs(x, a, b, metric = NULL)
```

## Arguments

- x:

  A `sim_results` object.

- a, b:

  `.config` labels.

- metric:

  Metric name; defaults to the first in the set.

## Value

A one-row tibble.
