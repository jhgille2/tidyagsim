# A soybean-style inbred line breeding pipeline

A ready-made
[`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md)
for the pipeline most public soybean (and other self-pollinated crop)
programs run: cross elite parents, advance the progeny to inbred lines
by single-seed descent, test every line in a preliminary yield trial
(PYT), advance the best to an advanced yield trial (AYT), release the
best AYT lines and recycle the best as parents for the next cycle. Any
argument can be tuned with
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html).

## Usage

``` r
line_pipeline(..., founders = asr_founders())
```

## Arguments

- ...:

  Arguments to set or tune: `n_parents`, `n_crosses`,
  `progeny_per_cross`, `line_generation`, `generations_per_year`,
  `pyt_locations`, `pyt_reps`, `n_ayt`, `ayt_locations`, `ayt_reps`,
  `n_release`, `cycles`, `gxe_var`, `cross_cost`, `ssd_cost`,
  `plot_cost`, `pyt_plot_cost`, `pyt_error`, `offseason_cost`,
  `max_plots`, `max_cost`, and any
  [`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md)
  field such as `plot_h2`.

- founders:

  An
  [`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md)
  recipe.

## Value

A `line_pipeline` spec.

## How one cycle works

1.  `n_crosses` random crosses among the parents (the founders in cycle
    1).

2.  Each F1 gives `progeny_per_cross` F2 plants, advanced by single-seed
    descent to F`line_generation` lines, so a cycle makes
    `n_crosses * progeny_per_cross` lines.

3.  PYT: every line at `pyt_locations` x `pyt_reps` plots. The top
    `n_ayt` lines by trial mean go on.

4.  AYT: `ayt_locations` x `ayt_reps` plots. The top `n_release` lines
    are released and the top `n_parents` become the next cycle's
    parents.

Trial means carry error variance
`gxe_var / locations + plot_var / (locations * reps)`, where `plot_var`
comes from the founders' `plot_h2` and `gxe_var` is in units of the
founder genetic variance, so extra locations average out
genotype-by-environment noise that extra reps cannot.

## Time and cost

A cycle takes `line_generation / generations_per_year` years in the
nursery (crossing and selfing) plus one year for each trial. Cycles do
not overlap, which keeps the model simple but understates the gain per
year of programs that start new crosses every year. Cost adds crosses
(`cross_cost`), single-seed-descent plants per generation (`ssd_cost`),
trial plots (`plot_cost`, or `pyt_plot_cost` for the PYT when it is
phenotyped differently) and, for nursery generations grown outside the
main season
(`line_generation - ceiling(line_generation / generations_per_year)` of
them), `offseason_cost` per line per generation. Designs whose plots or
cost per cycle exceed `max_plots` or `max_cost` are skipped as
infeasible.

## A different PYT phenotype

`pyt_error` multiplies the PYT error variance and `pyt_plot_cost` sets
its cost per plot, so a cheaper but noisier screen, such as drone
imagery or a visual score, can be compared with conventional yield
plots.

## Outcomes

`release_gain` (genetic value of the last cycle's releases, in founder
genetic SDs), `gain_per_year`, `gain_per_cost` (per 1000 currency
units), `genetic_variance` and `expected_heterozygosity` of the final
parents (relative to the founders), `accuracy` and `ayt_accuracy`
(correlation of PYT and AYT means with genetic value, among the lines in
each trial), `lines`, `plots` and `cost_per_cycle` per cycle, total
`cost`, `years`, and a yearly `trajectory` (see
[`collect_trajectories()`](https://jhgille2.github.io/tidyagsim/reference/collect_trajectories.md)).

## Examples

``` r
line_pipeline(n_ayt = tune(), ayt_locations = tune())
#> line_pipeline simulation spec (engine: AlphaSimR)
#> 
#> Tuned arguments:
#>   n_ayt = tune()
#>   ayt_locations = tune()
#> 
#> Fixed arguments:
#>   n_parents = 30
#>   n_crosses = 50
#>   progeny_per_cross = 20
#>   line_generation = 5
#>   generations_per_year = 2
#>   pyt_locations = 1
#>   pyt_reps = 1
#>   ayt_reps = 2
#>   n_release = 5
#>   cycles = 4
#>   gxe_var = 0.5
#>   cross_cost = 25
#>   ssd_cost = 2
#>   plot_cost = 15
#>   pyt_plot_cost = NA
#>   pyt_error = 1
#>   offseason_cost = 0
#>   max_plots = Inf
#>   max_cost = Inf
if (FALSE) { # \dontrun{
spec <- line_pipeline(n_ayt = tune(), ayt_locations = tune(), cycles = 3)
res <- tune_grid(spec, sim_replicates(10), grid = 9)
show_best(res, metric = "gain_per_year")
} # }
```
