# Pareto front of designs

Designs that no other design beats on every listed metric (comparing
means, respecting each metric's direction). The default front is the
first two metrics, for example power against cost.

## Usage

``` r
pareto_front(x, metrics = NULL, constraints = NULL, scenario = NULL)
```

## Arguments

- x:

  A `sim_results` object.

- metrics:

  Metric names; defaults to the first two in the set.

- constraints, scenario:

  As in
  [`show_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md).

## Value

A tibble with one row per design on the front, with means and standard
errors.
