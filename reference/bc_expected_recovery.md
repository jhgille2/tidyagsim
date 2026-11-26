# Backcrossing by the numbers

Classic expectations to check simulations against and to plan with.
`bc_expected_recovery()` gives the expected recurrent parent genome
share after `bc` backcross generations without background selection, 1 -
0.5^(bc + 1). `min_population_size()` gives the smallest number of
plants that contains at least one plant of a genotype with frequency
`p`, with probability `prob`: log(1 - prob) / log(1 - p).

## Usage

``` r
bc_expected_recovery(bc)

min_population_size(p, prob = 0.99)
```

## Arguments

- bc:

  Backcross generations (BC1 = 1).

- p:

  Frequency of the wanted genotype among progeny.

- prob:

  Wanted probability of finding at least one.

## Value

A numeric vector.

## Examples

``` r
bc_expected_recovery(1:6)
#> [1] 0.7500000 0.8750000 0.9375000 0.9687500 0.9843750 0.9921875
# Plants needed to find a BC1 carrier recombinant on both sides of a gene
# with flanking markers 10 cM away (Haldane):
r <- 0.5 * (1 - exp(-2 * 10 / 100))
min_population_size(0.5 * r^2)
#> [1] 1119
```
