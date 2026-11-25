# Run a tuner over every spec in a set

Like `workflowsets::workflow_map()`. Every spec uses the same
replicates, so designs from different specs are paired by seed.

## Usage

``` r
sim_map(object, fn = "tune_grid", ..., grids = NULL, verbose = FALSE)
```

## Arguments

- object:

  A
  [`sim_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_set.md).

- fn:

  Tuner name: `"tune_grid"`, `"tune_race"`, `"tune_bayes"` or
  `"fit_resamples"`.

- ...:

  Arguments passed to the tuner for every spec (e.g. `resamples`,
  `metrics`, `grid`).

- grids:

  Optional named list of grids, one per `wflow_id`, overriding `grid`.

- verbose:

  Print which spec is running.

## Value

The set with its `result` column filled.
