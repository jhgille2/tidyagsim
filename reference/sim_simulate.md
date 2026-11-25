# Run one design of a spec

The generic every simulator implements. It receives fully resolved
arguments and a vector of seeds and returns a list with one named list
of outcomes per seed (for example `list(power = 1, cost = 9800)`).

## Usage

``` r
sim_simulate(object, args, seeds, ...)
```

## Arguments

- object:

  A simulation spec.

- args:

  Resolved arguments (no
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) left).

- seeds:

  Integer seeds, one per replicate.

- ...:

  Not used.

## Value

A list of named lists, one per seed.
