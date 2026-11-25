# Collect tidy metrics

Like `tune::collect_metrics()`. With `summarize = TRUE`, one row per
design (and scenario, when there are several) and metric, with the mean
across replicates, the number of replicates `n` and the Monte Carlo
standard error `std_err`. When scenarios are pooled, each seed's value
is first averaged over scenarios so the standard error stays paired by
seed.

## Usage

``` r
collect_metrics(x, ...)

# S3 method for class 'sim_results'
collect_metrics(
  x,
  ...,
  summarize = TRUE,
  by_scenario = NULL,
  type = c("long", "wide")
)

# S3 method for class 'sim_set'
collect_metrics(
  x,
  ...,
  summarize = TRUE,
  by_scenario = NULL,
  type = c("long", "wide")
)
```

## Arguments

- x:

  A `sim_results` object.

- ...:

  Not used.

- summarize:

  Summarize over replicates (`TRUE`) or return one row per design,
  scenario and replicate (`FALSE`).

- by_scenario:

  Summarize within each scenario. Defaults to `TRUE` when the study has
  more than one scenario.

- type:

  `"long"` (one row per metric, like tune) or `"wide"` (one row per
  design, a column per metric and `<metric>_se` columns).

## Value

A tibble.
