# Mapping population design

A simulation spec for QTL mapping populations, so their design can be
tuned on its own or later as one stage of a breeding program. Each
replicate draws inbred founders, places one QTL of stated effect at a
random position, builds the population by simulated meiosis, phenotypes
the lines, scans the genome and records whether the QTL was detected.

## Usage

``` r
mapping_population(
  type = "RIL",
  lines = 200L,
  families = 10L,
  backcrosses = 1L,
  generation = 6L,
  founders = "4",
  locations = 1L,
  reps = 1L,
  qtl_pve = 0.05,
  plot_h2 = 0.5,
  gxe_ratio = 0.5,
  n_background = 1000L,
  qtl_donor_freq = 0.5,
  n_chromosomes = 20L,
  chromosome_length = 115,
  marker_spacing = 1,
  scan_model = "auto",
  n_perm = 100L,
  alpha = 0.05,
  threshold = NULL,
  window = 10,
  scan_step = NULL,
  cross_cost = 25,
  ssd_cost = 2,
  genotype_cost = 20,
  plot_cost = 15,
  generations_per_year = 2L,
  max_plots = Inf,
  max_cost = Inf,
  max_years = Inf
)
```

## Arguments

- type:

  Population type (see above).

- lines:

  Total lines phenotyped and genotyped.

- families:

  Number of families (NAM and AB-NAM).

- backcrosses:

  Backcross generations *s* (BCsFt and AB-NAM).

- generation:

  Filial generation *t* at phenotyping (RIL, BCsFt, NAM, AB-NAM, MAGIC).

- founders:

  MAGIC founders, `"4"` or `"8"`.

- locations, reps:

  Field locations and replicates per location.

- qtl_pve:

  Variance explained by the QTL in a reference RIL.

- plot_h2:

  Plot-basis heritability of the polygenic background.

- gxe_ratio:

  Genotype-by-location variance as a multiple of genetic variance.

- n_background:

  Number of background loci.

- qtl_donor_freq:

  Probability a donor or MAGIC founder carries the QTL allele.

- n_chromosomes, chromosome_length, marker_spacing:

  Genome (see
  [`sim_genetic_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_genetic_map.md)).

- scan_model:

  `"auto"`, `"marker"`, `"family_marker"` or `"haplotype"`.

- n_perm:

  Permutations per replicate for the threshold.

- alpha:

  Genome-wide type I error rate.

- threshold:

  Fixed LOD threshold (skips permutations; much faster).

- window:

  Distance in cM within which a peak counts as the QTL.

- scan_step:

  Scan every `scan_step` cM instead of every marker.

- cross_cost, ssd_cost, genotype_cost, plot_cost:

  Unit costs.

- generations_per_year:

  Generations per year (field, winter nursery, greenhouse).

- max_plots, max_cost, max_years:

  Budget limits; designs over them are infeasible.

## Value

A `mapping_population` simulation spec.

## Population types

Following R/qtl's BC*s*F*t* notation, with *t* the filial generation:

|  |  |  |
|----|----|----|
| `type` | Construction | Arguments used |
| `"F2"` | Two founders -\> F1 -\> one selfing (phenotyped as F2:3 family means) | `lines` |
| `"RIL"` | Two founders -\> F1 -\> single-seed descent to F*t* | `lines`, `generation` |
| `"BCsFt"` | Two founders -\> F1 -\> *s* backcrosses to founder 1 -\> self to BC*s*F*t* | `lines`, `backcrosses`, `generation` |
| `"NAM"` | One common founder x `families` donors -\> RILs in each family | `lines`, `families`, `generation` |
| `"AB-NAM"` | As NAM, with *s* backcrosses to the common parent in each family | `lines`, `families`, `backcrosses`, `generation` |
| `"MAGIC"` | 4 or 8 founders intermated in a funnel -\> self to F*t* | `lines`, `founders`, `generation` |

`lines` is always the total number of lines, so designs with the same
`lines`, `locations` and `reps` use the same field plots and compare
fairly. F1 is generation 1: F2 is `generation = 2`, an F6 RIL is
`generation = 6`, and BC1F3 is `backcrosses = 1, generation = 3` (one
backcross, then two selfings), matching the convention of breedsim's
`mapping-program`.

## QTL effect and genetic model

Traits are additive. `qtl_pve` is the share of single-plot phenotypic
variance the QTL would explain in a biparental RIL in which it
segregates, so its allele effect is the same in every design and
realized variance explained differs by design (it is reported as
`realized_pve`). In NAM and AB-NAM, donors carry the QTL allele with
probability `qtl_donor_freq`; in MAGIC every founder does. A polygenic
background of `n_background` loci gives the plot-basis heritability
`plot_h2`, and genotype-by-location variance is `gxe_ratio` times the
genetic variance, so more locations and reps sharpen line means.

## Detection

A replicate counts as a detection when the LOD within `window` cM of the
QTL exceeds the genome-wide threshold (permutations, or `threshold`).
`other_peaks` records whether any other chromosome also exceeds it. With
a polygenic background these are mostly real background loci (linked
small effects add up in segments that recombination has not broken), so
`other_peaks` rises with population size; with `n_background = 0` it is
a genome-wide false-positive rate, which the permutation threshold holds
near `alpha`. See
[`qtl_scan()`](https://jhgille2.github.io/tidyagsim/reference/qtl_scan.md)
for the models.

## Costs and time

`cost` adds crosses, single-seed-descent plant generations, genotyping
and field plots at the unit costs given (placeholders: replace them with
your own). `years` is generations divided by `generations_per_year`,
rounded up, plus one year for the field trial. `max_plots`, `max_cost`
and `max_years` make designs over budget infeasible, so tuners skip
them.

## Examples

``` r
spec <- mapping_population(type = "RIL", lines = tune(), generation = tune(), n_chromosomes = 4)
spec
#> mapping_population simulation spec (engine: tidyagsim)
#> 
#> Tuned arguments:
#>   lines = tune()
#>   generation = tune()
#> 
#> Fixed arguments:
#>   type = RIL
#>   families = 10
#>   backcrosses = 1
#>   founders = 4
#>   locations = 1
#>   reps = 1
#>   qtl_pve = 0.05
#>   plot_h2 = 0.5
#>   gxe_ratio = 0.5
#>   n_background = 1000
#>   qtl_donor_freq = 0.5
#>   n_chromosomes = 4
#>   chromosome_length = 115
#>   marker_spacing = 1
#>   scan_model = auto
#>   n_perm = 100
#>   alpha = 0.05
#>   threshold = NULL
#>   window = 10
#>   scan_step = NULL
#>   cross_cost = 25
#>   ssd_cost = 2
#>   genotype_cost = 20
#>   plot_cost = 15
#>   generations_per_year = 2
#>   max_plots = Inf
#>   max_cost = Inf
#>   max_years = Inf
extract_parameter_set_dials(spec)
#> Collection of 2 parameters for tuning
#> 
#>  identifier       type    object
#>       lines      lines nparam[+]
#>  generation generation nparam[+]
#> 
```
