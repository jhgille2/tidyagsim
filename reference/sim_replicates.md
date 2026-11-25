# Paired replicate seeds

In a simulation study the unit of replication is a random seed. Every
candidate design is run on the *same* seeds (common random numbers), so
differences between designs are paired: the same founders, QTL positions
and field errors are drawn for each design wherever the designs allow
it. `sim_replicates()` plays the role of `rsample::vfold_cv()`.

## Usage

``` r
sim_replicates(times = 10, seed = 1, scenarios = NULL)
```

## Arguments

- times:

  Number of replicates.

- seed:

  First seed; replicate `i` uses `seed + i - 1`.

- scenarios:

  Optional assumption scenarios from
  [`sim_scenarios()`](https://jhgille2.github.io/tidyagsim/reference/sim_scenarios.md).
  Each design is run on every replicate in every scenario.

## Value

A tibble of class `sim_rset` with columns `id` and `seed`, and the
scenarios as an attribute.

## Examples

``` r
sim_replicates(20)
#> # Simulation replicates: 20 paired seeds (1..20) x 1 scenario(s)
#> # A tibble: 20 × 2
#>    id     seed
#>    <chr> <int>
#>  1 Rep01     1
#>  2 Rep02     2
#>  3 Rep03     3
#>  4 Rep04     4
#>  5 Rep05     5
#>  6 Rep06     6
#>  7 Rep07     7
#>  8 Rep08     8
#>  9 Rep09     9
#> 10 Rep10    10
#> 11 Rep11    11
#> 12 Rep12    12
#> 13 Rep13    13
#> 14 Rep14    14
#> 15 Rep15    15
#> 16 Rep16    16
#> 17 Rep17    17
#> 18 Rep18    18
#> 19 Rep19    19
#> 20 Rep20    20
sim_replicates(10, scenarios = sim_scenarios(qtl_pve = c(0.03, 0.08)))
#> # Simulation replicates: 10 paired seeds (1..10) x 2 scenario(s)
#> # Scenarios: qtl_pve=0.03 | qtl_pve=0.08 
#> # A tibble: 10 × 2
#>    id     seed
#>    <chr> <int>
#>  1 Rep01     1
#>  2 Rep02     2
#>  3 Rep03     3
#>  4 Rep04     4
#>  5 Rep05     5
#>  6 Rep06     6
#>  7 Rep07     7
#>  8 Rep08     8
#>  9 Rep09     9
#> 10 Rep10    10
```
