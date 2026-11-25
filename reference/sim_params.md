# Parameter objects for simulation studies

dials parameter functions for arguments of
[`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md)
and other specs. Use them to change a tuned range, exactly as
`penalty()` is used in tidymodels:
`update(params, lines = n_lines(c(100, 600)))`.

## Usage

``` r
n_lines(range = c(50L, 1000L), trans = NULL)

n_families(range = c(2L, 25L), trans = NULL)

backcross_generations(range = c(0L, 3L), trans = NULL)

filial_generation(range = c(2L, 8L), trans = NULL)

population_type(values = c("F2", "RIL", "BCsFt", "NAM", "AB-NAM", "MAGIC"))

magic_founders(values = c("4", "8"))

n_locations(range = c(1L, 4L), trans = NULL)

n_reps(range = c(1L, 3L), trans = NULL)

qtl_pve(range = c(0.02, 0.15), trans = NULL)

plot_heritability(range = c(0.1, 0.8), trans = NULL)

generations_per_year(range = c(1L, 3L), trans = NULL)
```

## Arguments

- range:

  A two-element vector with the inclusive lower and upper bound.

- trans:

  A `trans` object from the scales package, or `NULL`.

- values:

  Allowed values of a qualitative parameter.

## Value

A dials `param` object.

## Examples

``` r
n_lines()
#> Lines in the population (quantitative)
#> Range: [50, 1000]
n_lines(c(100, 400))
#> Lines in the population (quantitative)
#> Range: [100, 400]
population_type(c("RIL", "NAM"))
#> Mapping population type (qualitative)
#> 2 possible values include:
#> 'RIL' and 'NAM'
```
