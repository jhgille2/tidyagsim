# Marker-assisted backcrossing and gene pyramiding

Simulates introgressing one or more target genes from a donor into a
recurrent parent (RP) by marker-assisted backcrossing, then fixing them
by selfing. Every argument can be tuned, so the classic design questions
(how many plants per generation, how many backcross generations, how
close the flanking markers should be, how dense the background markers,
and whether to introgress several genes together or in parallel streams)
become one call to
[`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md)
or
[`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md).

## Usage

``` r
backcross_program(..., targets = data.frame(chr = 1, pos = 50))
```

## Arguments

- ...:

  Arguments to set or tune: `strategy`, `bc_generations`, `n_plants`,
  `n_keep`, `flank_cM`, `bg_spacing`, `bg_from`, `fix_plants`,
  `pyramid_plants`, `n_chromosomes`, `chromosome_length`, `plant_cost`,
  `datapoint_cost`, `generations_per_year`.

- targets:

  Data frame with columns `chr` and `pos` (cM): where the target genes
  are.

## Value

A `backcross_program` spec.

## Selection in each backcross generation

Each generation grows `n_plants` plants and selects in three stages, as
in Hospital and Charcosset (1997) and Frisch et al. (1999):

1.  **Foreground:** plants carrying the donor allele at every target,
    read from a diagnostic marker on each gene.

2.  **Recombinant:** among carriers, those with the most recombinant
    flanking markers, placed `flank_cM` on each side of each target.
    This limits linkage drag. Use `flank_cM = Inf` to skip this stage.

3.  **Background:** from generation `bg_from` on, the remaining plants
    are ranked by the share of RP alleles at markers every `bg_spacing`
    cM. Use `bg_spacing = Inf` to skip it.

The best `n_keep` plants are backcrossed to the RP. After
`bc_generations`, the best plant is selfed, `fix_plants` progeny are
grown, and a plant homozygous for every target is kept, with ties broken
by background.

## Several targets

With more than one target, `strategy = "simultaneous"` backcrosses all
of them together, selecting carriers of every target each generation.
`strategy = "parallel"` follows Peng et al. (2014a, b). Each target gets
its own backcross stream of `n_plants` per generation, and each
conversion is fixed by selfing. The conversions are then pyramided by
crossing them in pairs, then crossing those hybrids, and so on (a
symmetric schedule, as in Ishii and Yonezawa, 2007). Each crossing step
grows `pyramid_plants` and keeps the plant carrying every target with
the best background. A final selfing fixes the pyramid.

## Outcomes

- `success`: 1 if the final plant is homozygous for the donor allele at
  every target, else 0.

- `rp_recovery`: percent of the final plant's genome from the RP.

- `donor_cM`: donor genome left, in cM per haploid genome.

- `linkage_drag`: mean length in cM of the donor segment around each
  target.

- `generations` and `years`: crossing and selfing generations, and years
  at `generations_per_year`. Parallel streams run at the same time.

- `plants`, `datapoints` (marker data points) and `cost`.

## References

Hospital F, Charcosset A (1997) Marker-assisted introgression of
quantitative trait loci. Genetics 147:1469-1485.
[doi:10.1093/genetics/147.3.1469](https://doi.org/10.1093/genetics/147.3.1469)

Frisch M, Bohn M, Melchinger AE (1999) Comparison of selection
strategies for marker-assisted backcrossing of a gene. Crop Science
39:1295-1301.

Hospital F (2001) Size of donor chromosome segments around introgressed
loci and reduction of linkage drag in marker-assisted backcross
programs. Genetics 158:1363-1379.
[doi:10.1093/genetics/158.3.1363](https://doi.org/10.1093/genetics/158.3.1363)

Peng T, Sun X, Mumm RH (2014a) Optimized breeding strategies for
multiple trait integration: I. Minimizing linkage drag in single event
introgression. Molecular Breeding 33:89-104.
[doi:10.1007/s11032-013-9936-7](https://doi.org/10.1007/s11032-013-9936-7)

Peng T, Sun X, Mumm RH (2014b) Optimized breeding strategies for
multiple trait integration: II. Process efficiency in event pyramiding
and trait fixation. Molecular Breeding 33:105-115.
[doi:10.1007/s11032-013-9937-6](https://doi.org/10.1007/s11032-013-9937-6)

Ishii T, Yonezawa K (2007) Optimization of the marker-based procedures
for pyramiding genes from multiple donor lines: I. Schedule of crossing
between the donor lines. Crop Science 47:537-546.

## Examples

``` r
spec <- backcross_program(n_plants = tune(), flank_cM = tune(), n_chromosomes = 5)
res <- tune_grid(spec, sim_replicates(5), grid = expand.grid(n_plants = c(50, 150), flank_cM = c(5, 20)))
collect_metrics(res, type = "wide")
#> # A tibble: 4 × 18
#>   n_plants flank_cM success success_se rp_recovery rp_recovery_se linkage_drag
#>      <dbl>    <dbl>   <dbl>      <dbl>       <dbl>          <dbl>        <dbl>
#> 1       50        5       1          0        98.8          0.191          3.4
#> 2      150        5       1          0        98.6          0.265          5  
#> 3       50       20       1          0        98.0          0.463          9.8
#> 4      150       20       1          0        97.7          0.340         12.4
#> # ℹ 11 more variables: linkage_drag_se <dbl>, years <dbl>, years_se <dbl>,
#> #   plants <dbl>, plants_se <dbl>, datapoints <dbl>, datapoints_se <dbl>,
#> #   cost <dbl>, cost_se <dbl>, n <int>, .config <chr>
```
