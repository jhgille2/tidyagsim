# Race designs and drop the clearly worse ones early

Like `finetune::tune_race_anova()`. Every design runs on the first
`burn_in` replicates. After each further replicate, a design is dropped
when a one-sided paired t-test on shared seeds says it is worse than the
current leader at level `alpha`. Designs that break `constraints` on
current means are dropped as soon as a design meeting them leads.
Survivors are run on every replicate, so their estimates are comparable
with a full grid.

## Usage

``` r
tune_race(
  object,
  resamples = NULL,
  ...,
  grid = NULL,
  metrics = NULL,
  param_info = NULL,
  control = control_sim(),
  burn_in = 3,
  alpha = 0.05,
  metric = NULL,
  constraints = NULL,
  seed = 1
)
```

## Arguments

- object:

  A simulation spec with
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html)
  arguments.

- resamples:

  Replicates from
  [`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md),
  or a number of replicates.

- ...:

  Not used.

- grid:

  `NULL` (10 space-filling designs), a number of designs, or a data
  frame of designs with one column per tuned parameter id.

- metrics:

  A
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md);
  defaults to the spec's own.

- param_info:

  A dials `parameters` object with changed ranges, from
  `extract_parameter_set_dials(spec) |> update(...)`.

- control:

  A
  [`control_sim()`](https://jhgille2.github.io/tidyagsim/reference/control_sim.md)
  object.

- burn_in:

  Replicates every design gets before any is dropped (at least 2).

- alpha:

  Significance level for dropping a design.

- metric:

  The metric to race on; defaults to the first in the set.

- constraints:

  Named list of bounds on other metrics' means, e.g.
  `list(power = 0.8)`. A bare number is a floor for metrics to maximize
  and a ceiling for metrics to minimize; `list(cost = c("<=", 5000))` is
  explicit.

- seed:

  Seed for generating the space-filling grid.

## Value

A `sim_results` object; `x$trace` records how many designs survived each
round.
