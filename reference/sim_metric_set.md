# Combine metrics into a metric set

Like `yardstick::metric_set()`. Name built-in metrics bare or as strings
and add custom ones made with
[`sim_metric()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric.md).
The first metric is the one tuners optimize and selection helpers rank
by unless told otherwise.

## Usage

``` r
sim_metric_set(...)
```

## Arguments

- ...:

  Metric names (bare or quoted) or `sim_metric` objects.

## Value

A `sim_metric_set`.

## Details

Built-in metrics: `power`, `other_peaks`, `location_error`,
`support_interval`, `peak_lod`, `realized_pve`, `threshold` (mapping);
`release_gain`, `gain_per_year`, `gain_per_cost`,
`expected_heterozygosity`, `genetic_variance` (breeding programs);
`lines`, `plots`, `genotyped`, `crosses`, `cost`, `years`, `generations`
(resources).

## Examples

``` r
sim_metric_set(power, cost, years)
#> A simulation metric set:
#>   power              maximize  QTL detection power
#>   cost               minimize  Total cost
#>   years              minimize  Years from first cross to QTL results
```
