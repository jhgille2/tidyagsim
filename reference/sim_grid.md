# Build the candidate grid for a study

`NULL` or a number gives that many space-filling designs (10 by
default), or every combination when the space is that small. A data
frame is used as given after checking its columns. `expand_grid`-style
tibbles from
[`dials::grid_regular()`](https://dials.tidymodels.org/reference/grid_regular.html)
work directly.

## Usage

``` r
sim_grid(params, grid = NULL, seed = 1)
```

## Arguments

- params:

  A dials `parameters` object.

- grid:

  `NULL`, a number of designs, or a data frame of designs.

- seed:

  Seed for the space-filling design.

## Value

A tibble with `.config` and one column per parameter id.
