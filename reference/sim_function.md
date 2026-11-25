# A simulation spec from any R function

Wraps a function `fn(..., seed)` that runs one replicate and returns a
named list (or one-row data frame) of outcomes, or a single number
(stored as `value`). Arguments given here are passed to `fn`; set any of
them to [`tune()`](https://hardhat.tidymodels.org/reference/tune.html)
and give its range in `params`.

## Usage

``` r
sim_function(fn, ..., params = list(), check = NULL, metrics = NULL)
```

## Arguments

- fn:

  A function whose arguments include the ones given in `...` and `seed`.

- ...:

  Arguments for `fn`; use
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) to mark
  one for tuning.

- params:

  Named list of dials parameter objects giving default ranges for tuned
  arguments (see
  [`sim_param()`](https://jhgille2.github.io/tidyagsim/reference/sim_param.md)).

- check:

  Optional function of the same arguments returning `NULL` when a design
  is feasible or a reason string when it is not.

- metrics:

  Optional default
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md).

## Value

A `sim_function` spec.

## Examples

``` r
toy <- sim_function(
  function(x, noise, seed) {
    set.seed(seed)
    list(y = -(x - 3)^2 + rnorm(1, sd = noise))
  },
  x = tune(), noise = 0.5,
  params = list(x = sim_param(c(0, 6)))
)
res <- tune_grid(toy, sim_replicates(5), grid = 7,
                 metrics = sim_metric_set(sim_metric("y", "maximize")))
show_best(res, n = 2)
#> # A tibble: 2 × 6
#>       x .metric   mean     n std_err .config  
#>   <int> <chr>    <dbl> <int>   <dbl> <chr>    
#> 1     3 y       -0.311     5   0.109 Config004
#> 2     2 y       -1.31      5   0.109 Config003
```
