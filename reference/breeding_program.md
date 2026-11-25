# A breedsim breeding program as a simulation spec

Wraps a flowchart program drawn in the breedsim desktop editor (or
written as JSON) so its fields can be tuned with the tidyagsim verbs.
Arguments are named by `"stage.field"` paths, e.g.
`` `pyt.reps` = tune() ``. Fields written as `{"tune": ...}` in the
program JSON are marked for tuning automatically, with their bounds.
Each replicate is simulated by the AlphaSimR engine that ships with
breedsim.

## Usage

``` r
breeding_program(
  path,
  ...,
  year = NULL,
  start_year = NULL,
  parent_stage = NULL,
  discount_rate = 0,
  settings = NULL,
  rscript = NULL
)
```

## Arguments

- path:

  A breedsim program or study JSON file.

- ...:

  `"stage.field" = value` or
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html).

- year:

  Release year scored (defaults to the program horizon).

- start_year:

  First year of the window for rate metrics.

- parent_stage:

  A crossing stage whose parent-pool diversity is reported.

- discount_rate:

  Annual discount rate for `gain_per_cost`.

- settings:

  Path to a breedsim simulation settings JSON file.

- rscript:

  Path to `Rscript` for the engine (defaults to the one on `PATH` or the
  `BREEDSIM_RSCRIPT` environment variable).

## Value

A `breeding_program` simulation spec.

## Details

Reported outcomes per replicate: `release_gain` in `year`;
`gain_per_year` and `gain_per_cost` over `start_year`..`year` when
`start_year` is given; `expected_heterozygosity` of the crossing pool at
`parent_stage` when given; and the program's scheduled `plots`,
`genotyped` and `cost` summed over the planning horizon.

## Examples

``` r
if (FALSE) { # \dontrun{
spec <- breeding_program("examples/conventional.json",
                         `pyt.reps` = tune(), `pyt.locations` = tune())
res <- tune_grid(spec, sim_replicates(5), grid = 8,
                 metrics = sim_metric_set(release_gain, plots))
} # }
```
