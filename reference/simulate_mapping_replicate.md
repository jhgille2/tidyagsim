# Simulate one mapping population replicate

The engine behind
[`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md).
Exposed for checking and teaching: it returns the genotypes, phenotypes,
QTL and scan as well as the outcomes.

## Usage

``` r
simulate_mapping_replicate(args, seed, keep_data = FALSE)
```

## Arguments

- args:

  A full argument list, as built by
  [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md).

- seed:

  Replicate seed.

- keep_data:

  Return the simulated genotypes, phenotypes and LOD profile.

## Value

A list of outcomes (`power` is 0 or 1 for one replicate), plus data when
`keep_data = TRUE`.
