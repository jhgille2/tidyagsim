# Define a simulation metric

A metric reads one number from a replicate's outcomes and says whether
bigger or smaller is better. Built-in metrics (see
[`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md))
need no definition; use `sim_metric()` for anything else.

## Usage

``` r
sim_metric(name, direction = "maximize", fn = NULL, label = name, units = "")
```

## Arguments

- name:

  Metric name; also the outcome to read when `fn` is `NULL`.

- direction:

  `"maximize"`, `"minimize"` or `"zero"` (closest to zero is best).

- fn:

  Optional function of the outcome list returning one number, for a
  derived metric.

- label, units:

  Labels used in printing and plots.

## Value

A `sim_metric` object.

## Examples

``` r
cost_per_line <- sim_metric("cost_per_line", "minimize", fn = function(o) o$cost / o$lines)
sim_metric_set(power, cost, cost_per_line)
#> A simulation metric set:
#>   power              maximize  QTL detection power
#>   cost               minimize  Total cost
#>   cost_per_line      minimize  cost_per_line
```
