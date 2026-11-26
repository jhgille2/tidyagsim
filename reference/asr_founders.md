# Founder population recipe for AlphaSimR

Describes how each replicate builds its founders and trait with
[AlphaSimR](https://cran.r-project.org/package=AlphaSimR). Every
replicate draws new founders from its own seed, and every design in a
study uses the same seeds, so all designs start from the same founders
(common random numbers).

## Usage

``` r
asr_founders(
  n_founders = 50,
  n_chr = 20,
  seg_sites = 200,
  qtl_per_chr = 20,
  mean = 0,
  var = 1,
  plot_h2 = 0.3,
  species = "none",
  haplotypes = NULL,
  setup = NULL
)
```

## Arguments

- n_founders:

  Number of inbred founder lines.

- n_chr:

  Number of chromosomes (soybean has 20).

- seg_sites:

  Segregating sites per chromosome.

- qtl_per_chr:

  QTL per chromosome for the additive trait.

- mean, var:

  Founder mean and genetic variance of the trait.

- plot_h2:

  Heritability of a single plot in the founders; sets the plot error
  variance.

- species:

  `"none"` or a MaCS species history; see below.

- haplotypes:

  Optional AlphaSimR `MapPop` to sample founders from.

- setup:

  Optional function of `SP` (the `SimParam`) for anything else, such as
  extra traits (`SP$addTraitAD()`), genetic correlations or SNP chips.
  When given, it replaces the default trait.

## Value

An `asr_founders` recipe.

## Details

Any of these fields except `setup` and `haplotypes` can also be given to
[`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md)
or
[`line_pipeline()`](https://jhgille2.github.io/tidyagsim/reference/line_pipeline.md)
as an argument, so they can be tuned or turned into scenarios, for
example `sim_scenarios(plot_h2 = c(0.1, 0.3))`.

## Where founder haplotypes come from

- `species = "none"` (the default) draws founder haplotypes with
  [`AlphaSimR::quickHaplo()`](https://gaynorr.github.io/AlphaSimR/reference/quickHaplo.html):
  fast and exactly reproducible from the seed, but without the linkage
  disequilibrium a real population history leaves.

- `species = "GENERIC"`, `"WHEAT"`, `"MAIZE"`, ... simulates a
  population history with MaCS
  ([`AlphaSimR::runMacs()`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.html)).
  MaCS does not follow R's random seed, so tidyagsim keeps each
  replicate's founders in memory for the session: designs in one session
  share founders, but a new session (or a parallel worker) draws
  different ones.

- `haplotypes` takes a `MapPop` you made yourself, for example from real
  marker data with
  [`AlphaSimR::importHaplo()`](https://gaynorr.github.io/AlphaSimR/reference/importHaplo.html)
  or one MaCS run saved with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Each replicate
  samples `n_founders` of its individuals.

## Examples

``` r
asr_founders(n_founders = 40, plot_h2 = 0.2)
#> # AlphaSimR founders: 40 inbred lines, 20 chromosomes x 200 sites (independent haplotypes)
#> # Trait: 20 QTL per chromosome, plot h2 0.2
```
