# Assumption scenarios

Crosses assumption factors into named scenarios, for things a breeder
cannot choose but must plan around (heritability, marker cost, QTL
effect size). Each name must be an argument of the spec being tuned; the
scenario value overrides the spec's fixed value.

## Usage

``` r
sim_scenarios(...)
```

## Arguments

- ...:

  Named vectors. Scalars are held fixed but still recorded.

## Value

A tibble of class `sim_scenarios` with a `.scenario` label column and
one column per factor.

## Examples

``` r
sim_scenarios(plot_h2 = c(0.3, 0.6), genotype_cost = c(15, 30))
#> # A tibble: 4 × 3
#>   .scenario                     plot_h2 genotype_cost
#> * <chr>                           <dbl>         <dbl>
#> 1 plot_h2=0.3, genotype_cost=15     0.3            15
#> 2 plot_h2=0.6, genotype_cost=15     0.6            15
#> 3 plot_h2=0.3, genotype_cost=30     0.3            30
#> 4 plot_h2=0.6, genotype_cost=30     0.6            30
```
