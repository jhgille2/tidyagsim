# Simulation specifications

A simulation spec is the tidyagsim analogue of a parsnip model
specification. It records the arguments of one simulator. Any argument
may be set to
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html) to let a
tuner choose it; the rest stay fixed. Every spec can be passed to
[`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md),
[`tune_race()`](https://jhgille2.github.io/tidyagsim/reference/tune_race.md),
[`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md)
and
[`fit_resamples()`](https://jhgille2.github.io/tidyagsim/reference/fit_resamples.md).

## Usage

``` r
new_sim_spec(cls, args, engine = "tidyagsim", ...)
```

## Arguments

- cls:

  Subclass name.

- args:

  Named list of arguments; values may be
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) calls.

- engine:

  Engine name (the simulator implementation).

- ...:

  Other fields stored on the spec.

## Value

An object of class `c(cls, "sim_spec")`.

## Details

Spec constructors shipped with the package:

- [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md):
  QTL mapping populations, simulated in R.

- [`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md):
  any R function of the arguments and a seed.

- [`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md):
  a breedsim flowchart program run by AlphaSimR.

New simulators are added by writing a constructor that calls
`new_sim_spec()` and methods for
[`sim_simulate()`](https://jhgille2.github.io/tidyagsim/reference/sim_simulate.md)
and, optionally,
[`sim_check()`](https://jhgille2.github.io/tidyagsim/reference/sim_check.md)
and
[`sim_parameters()`](https://jhgille2.github.io/tidyagsim/reference/sim_parameters.md).
