# Collect every outcome the simulator reported

Needs `control_sim(save_outcomes = TRUE)`. Returns one row per design,
scenario and replicate, with one column per scalar outcome (for example
the QTL position or the number of segregating markers in a mapping
population), like `tune::collect_predictions()`.

## Usage

``` r
collect_outcomes(x)
```

## Arguments

- x:

  A `sim_results` object.

## Value

A tibble.
