# Run every design in a grid

The simulation analogue of `tune::tune_grid()`: every design runs on
every replicate seed and scenario.

## Usage

``` r
tune_grid(object, ...)

# S3 method for class 'sim_spec'
tune_grid(
  object,
  resamples = NULL,
  ...,
  grid = NULL,
  metrics = NULL,
  param_info = NULL,
  control = control_sim(),
  seed = 1
)
```

## Arguments

- object:

  A simulation spec with
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html)
  arguments.

- ...:

  Not used.

- resamples:

  Replicates from
  [`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md),
  or a number of replicates.

- grid:

  `NULL` (10 space-filling designs), a number of designs, or a data
  frame of designs with one column per tuned parameter id.

- metrics:

  A
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md);
  defaults to the spec's own.

- param_info:

  A dials `parameters` object with changed ranges, from
  `extract_parameter_set_dials(spec) |> update(...)`.

- control:

  A
  [`control_sim()`](https://jhgille2.github.io/tidyagsim/reference/control_sim.md)
  object.

- seed:

  Seed for generating the space-filling grid.

## Value

A `sim_results` object; see
[`collect_metrics()`](https://jhgille2.github.io/tidyagsim/reference/collect_metrics.md).

## Examples

``` r
spec <- mapping_population(type = "RIL", lines = tune(), n_chromosomes = 2, n_perm = 20)
res <- tune_grid(spec, sim_replicates(3), grid = data.frame(lines = c(60, 120)),
                 metrics = sim_metric_set(power, cost))
collect_metrics(res)
#> # A tibble: 4 × 7
#>   lines .metric direction     mean     n std_err .config  
#>   <dbl> <chr>   <chr>        <dbl> <int>   <dbl> <chr>    
#> 1    60 power   maximize     0         3   0     Config001
#> 2    60 cost    minimize  2725         3   0     Config001
#> 3   120 power   maximize     0.333     3   0.333 Config002
#> 4   120 cost    minimize  5425         3   0     Config002
```
