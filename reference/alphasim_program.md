# A simulation spec from an AlphaSimR script

Turns an AlphaSimR breeding scheme into a tunable spec. You write one
replicate of the scheme as a function of the founders, the `SimParam`
and your design arguments; tidyagsim builds the founders for each
replicate, seeds the random numbers, runs every design on the same seeds
and collects the results. Mark any argument with
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html) to search
over it.

## Usage

``` r
alphasim_program(
  program,
  ...,
  founders = asr_founders(),
  params = list(),
  check = NULL,
  metrics = NULL
)
```

## Arguments

- program:

  A function whose first two arguments are `founders` (an AlphaSimR
  `Pop`) and `SP`, followed by design arguments. It returns a named list
  of outcomes.

- ...:

  Values for `program`'s arguments and for any
  [`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md)
  field; use
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) to mark
  one for tuning.

- founders:

  An
  [`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md)
  recipe.

- params:

  Named list of parameter objects (see
  [`sim_param()`](https://jhgille2.github.io/tidyagsim/reference/sim_param.md))
  giving ranges for tuned arguments.

- check:

  Optional function of the design arguments returning `NULL` when a
  design is feasible or a reason string when it is not.

- metrics:

  Optional default
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md).

## Value

An `alphasim_program` spec.

## Details

Inside `program`, `SP` is also bound globally while it runs, so
AlphaSimR calls work with or without `simParam = SP`.

Use
[`asr_tracker()`](https://jhgille2.github.io/tidyagsim/reference/asr_tracker.md)
inside `program` to record the population each year and
[`asr_summary()`](https://jhgille2.github.io/tidyagsim/reference/asr_summary.md)
to turn the record into the built-in breeding metrics (`release_gain`,
`gain_per_year`, `genetic_variance`, `cost`, ...).

## Examples

``` r
if (FALSE) { # \dontrun{
recurrent <- alphasim_program(
  function(founders, SP, n_select, cycles) {
    tr <- asr_tracker(founders)
    pop <- founders
    for (cycle in seq_len(cycles)) {
      pop <- AlphaSimR::setPheno(pop)
      parents <- AlphaSimR::selectInd(pop, n_select)
      pop <- AlphaSimR::randCross(parents, nCrosses = 100)
      tr$record(pop, year = cycle, cost = 100 * 15)
    }
    asr_summary(tr)
  },
  n_select = tune(), cycles = 5,
  founders = asr_founders(n_chr = 5, seg_sites = 100),
  params = list(n_select = sim_param(c(5, 40)))
)
res <- tune_grid(recurrent, sim_replicates(10), grid = 5)
autoplot(res)
} # }
```
