# Plot simulation results

- `type = "marginals"`: mean +/- one standard error of each metric
  against each numeric parameter (qualitative parameters are shown as
  colour).

- `type = "pareto"`: the first two metrics against each other, with the
  Pareto front highlighted.

- `type = "race"` (racing results): designs remaining after each
  replicate.

## Usage

``` r
# S3 method for class 'sim_results'
autoplot(object, type = c("marginals", "pareto", "race"), metric = NULL, ...)

# S3 method for class 'sim_set'
autoplot(object, type = c("pareto", "ranks"), metric = NULL, ...)
```

## Arguments

- object:

  A `sim_results` object or a run
  [`sim_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_set.md).

- type:

  Plot type.

- metric:

  Metric names to show; defaults to all (marginals) or the first two
  (pareto).

- ...:

  Not used.

## Value

A ggplot.
