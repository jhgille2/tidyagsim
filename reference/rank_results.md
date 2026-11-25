# Rank designs across a set of specs

Like `workflowsets::rank_results()`.

## Usage

``` r
rank_results(x, rank_metric = NULL, constraints = NULL, select_best = FALSE)
```

## Arguments

- x:

  A
  [`sim_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_set.md)
  that has been run with
  [`sim_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_map.md).

- rank_metric:

  Metric to rank by; defaults to the first in the set.

- constraints:

  Named list of bounds on metric means, as in
  [`show_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md).

- select_best:

  Keep only the best design of each spec.

## Value

A wide tibble, best first, with a `rank` column.
