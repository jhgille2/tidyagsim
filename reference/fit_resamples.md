# Run one fixed design on every replicate

The simulation analogue of `tune::fit_resamples()`.

## Usage

``` r
fit_resamples(object, ...)

# S3 method for class 'sim_spec'
fit_resamples(
  object,
  resamples = NULL,
  ...,
  metrics = NULL,
  control = control_sim()
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

- metrics:

  A
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md);
  defaults to the spec's own.

- control:

  A
  [`control_sim()`](https://jhgille2.github.io/tidyagsim/reference/control_sim.md)
  object.

## Value

A `sim_results` object with one design.

## Examples

``` r
spec <- mapping_population(type = "F2", lines = 100, n_chromosomes = 2, n_perm = 20)
collect_metrics(fit_resamples(spec, 3, metrics = sim_metric_set(power)))
#> # A tibble: 1 × 6
#>   .metric direction  mean     n std_err .config  
#>   <chr>   <chr>     <dbl> <int>   <dbl> <chr>    
#> 1 power   maximize      0     3       0 Config001
```
