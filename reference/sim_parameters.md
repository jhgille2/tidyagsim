# Default parameter objects of a spec

Returns a named list of dials parameter objects, one per argument that
has a sensible default range. Spec authors implement this generic so
that [`tune()`](https://hardhat.tidymodels.org/reference/tune.html)
works without the user writing ranges.

## Usage

``` r
sim_parameters(object, ...)
```

## Arguments

- object:

  A simulation spec.

- ...:

  Not used.

## Value

A named list of `param` objects.
