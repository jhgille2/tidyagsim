# Record a breeding program as it runs

Create a tracker from the founders at the start of an
[`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md),
call `$record()` with the population of interest each year (and the
money spent that year), then pass the tracker to
[`asr_summary()`](https://jhgille2.github.io/tidyagsim/reference/asr_summary.md).
Genetic values are reported in founder genetic standard deviations, so
gains are comparable across traits and scenarios.

## Usage

``` r
asr_tracker(founders, trait = 1)
```

## Arguments

- founders:

  The founder `Pop`.

- trait:

  Trait number to follow.

## Value

An `asr_tracker`: a list of functions
`record(pop, year, cost = 0, released_lines = NULL)`, `add_cost(cost)`
and `trajectory()`.
