# Collect yearly trajectories from a study

For specs whose outcomes include a `trajectory` data frame, such as
[`line_pipeline()`](https://jhgille2.github.io/tidyagsim/reference/line_pipeline.md)
or an
[`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md)
that returns
[`asr_summary()`](https://jhgille2.github.io/tidyagsim/reference/asr_summary.md).
Needs `control_sim(save_outcomes = TRUE)`.

## Usage

``` r
collect_trajectories(x)
```

## Arguments

- x:

  A `sim_results` object.

## Value

A tibble with one row per design, scenario, replicate and year,
including the tuned parameter values.
