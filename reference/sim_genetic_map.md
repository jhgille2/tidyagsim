# Simulate a genetic map

Evenly spaced markers on chromosomes of equal genetic length. The
defaults (20 chromosomes of 115 cM, a marker per cM) approximate the
size of published soybean consensus linkage maps and a mid-density SNP
panel.

## Usage

``` r
sim_genetic_map(
  n_chromosomes = 20,
  chromosome_length = 115,
  marker_spacing = 1
)
```

## Arguments

- n_chromosomes:

  Number of chromosomes.

- chromosome_length:

  Genetic length of each chromosome in cM (recycled).

- marker_spacing:

  Distance between markers in cM.

## Value

A list with `chr`, `pos` (cM) and `n`.

## Examples

``` r
map <- sim_genetic_map(2, 50, 5)
str(map)
#> List of 4
#>  $ chr    : int [1:22] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ pos    : num [1:22] 0 5 10 15 20 25 30 35 40 45 ...
#>  $ n      : int 22
#>  $ lengths: num [1:2] 50 50
```
