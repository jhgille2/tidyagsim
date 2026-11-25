# Find the best designs

`show_best()` lists the top designs by the mean of one metric;
`select_best()` returns the parameter values of the best one;
`select_by_one_std_err()` returns the simplest design within one
standard error of the best, and `select_by_pct_loss()` the simplest
within a percentage loss. These mirror the tune functions of the same
names.

## Usage

``` r
show_best(x, ...)

# S3 method for class 'sim_results'
show_best(x, ..., metric = NULL, n = 5, constraints = NULL, scenario = NULL)

select_best(x, ...)

# S3 method for class 'sim_results'
select_best(x, ..., metric = NULL, constraints = NULL, scenario = NULL)

select_by_one_std_err(x, ...)

# S3 method for class 'sim_results'
select_by_one_std_err(
  x,
  ...,
  metric = NULL,
  constraints = NULL,
  scenario = NULL
)

select_by_pct_loss(x, ...)

# S3 method for class 'sim_results'
select_by_pct_loss(
  x,
  ...,
  metric = NULL,
  limit = 2,
  constraints = NULL,
  scenario = NULL
)
```

## Arguments

- x:

  A `sim_results` object.

- ...:

  For `select_by_one_std_err()` and `select_by_pct_loss()`: columns
  (parameters or metrics) that order designs from simplest, e.g. `lines`
  or `desc(families)` or `cost`.

- metric:

  Metric to rank by; defaults to the first in the set.

- n:

  Number of designs to show.

- constraints:

  Named list of bounds on metric means.

- scenario:

  A scenario label to rank within, or `NULL` to pool.

- limit:

  Allowed percentage loss for `select_by_pct_loss()`.

## Value

A tibble.

## Details

Means are pooled over scenarios unless `scenario` names one. Constraints
restrict which designs are eligible, e.g.
`constraints = list(power = 0.8)` (at least 0.8, because power is
maximized) or `list(cost = c("<=", 20000))`.
