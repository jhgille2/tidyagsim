# Change the arguments of a spec

Change the arguments of a spec

## Usage

``` r
set_args(object, ...)
```

## Arguments

- object:

  A simulation spec.

- ...:

  Named arguments; use
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) to mark
  one for tuning.

## Value

The modified spec.

## Examples

``` r
mapping_population(type = "F2") |> set_args(lines = tune())
#> mapping_population simulation spec (engine: tidyagsim)
#> 
#> Tuned arguments:
#>   lines = tune()
#> 
#> Fixed arguments:
#>   type = F2
#>   families = 10
#>   backcrosses = 1
#>   generation = 6
#>   founders = 4
#>   locations = 1
#>   reps = 1
#>   qtl_pve = 0.05
#>   plot_h2 = 0.5
#>   gxe_ratio = 0.5
#>   n_background = 1000
#>   qtl_donor_freq = 0.5
#>   n_chromosomes = 20
#>   chromosome_length = 115
#>   marker_spacing = 1
#>   scan_model = auto
#>   n_perm = 100
#>   alpha = 0.05
#>   threshold = NULL
#>   window = 10
#>   scan_step = NULL
#>   cross_cost = 25
#>   ssd_cost = 2
#>   genotype_cost = 20
#>   plot_cost = 15
#>   generations_per_year = 2
#>   max_plots = Inf
#>   max_cost = Inf
#>   max_years = Inf
```
