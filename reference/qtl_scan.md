# Genome scan for a single QTL

Single-marker regression (`"marker"`), marker-within-family regression
as in nested association mapping joint linkage (`"family_marker"`), or
regression on founder-origin dosages (`"haplotype"`, for MAGIC). LOD
scores compare the residual sum of squares with and without the locus.
The genome-wide threshold is the `1 - alpha` quantile of the maximum LOD
over `n_perm` permutations of the phenotypes (within families when there
are several).

## Usage

``` r
qtl_scan(
  markers,
  y,
  family,
  map,
  model = c("marker", "family_marker", "haplotype"),
  n_perm = 100,
  alpha = 0.05,
  threshold = NULL,
  step = NULL,
  origin = NULL
)
```

## Arguments

- markers:

  Marker dosage matrix (markers x lines, coded 0/1/2).

- y:

  Phenotypes.

- family:

  Family of each line.

- map:

  Genetic map from
  [`sim_genetic_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_genetic_map.md).

- model:

  Scan model.

- n_perm:

  Permutations for the threshold.

- alpha:

  Genome-wide type I error rate.

- threshold:

  A fixed LOD threshold instead of permutations.

- step:

  Scan every `step` cM instead of every marker (`NULL` for every
  marker).

- origin:

  Founder origin matrices for the haplotype model.

## Value

A list with the scan positions, `lod` and `threshold`.

## Details

The haplotype model uses the true founder origin of each allele, so it
is an upper bound on what haplotype reconstruction (e.g. with qtl2)
achieves.
