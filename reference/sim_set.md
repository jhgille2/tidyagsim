# A set of simulation specs to compare

Like `workflowsets::workflow_set()`. Useful when candidates take
different parameters, for example RIL designs (tune `selfing`) versus
NAM designs (tune `families`) versus MAGIC (tune `founders`): each spec
gets its own grid, but all run on the same replicate seeds and are
ranked together.

## Usage

``` r
sim_set(...)
```

## Arguments

- ...:

  Named simulation specs.

## Value

A tibble of class `sim_set` with columns `wflow_id`, `spec` and
`result`.

## Examples

``` r
sim_set(
  ril = mapping_population(type = "RIL", lines = tune()),
  nam = mapping_population(type = "NAM", lines = tune(), families = tune())
)
#> # A simulation set of 2 spec(s)
#> # A tibble: 2 × 4
#>   wflow_id spec               tuned           result   
#>   <chr>    <chr>              <chr>           <chr>    
#> 1 ril      mapping_population lines           <not run>
#> 2 nam      mapping_population lines, families <not run>
```
