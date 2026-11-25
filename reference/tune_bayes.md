# Bayesian optimization of a simulation study

Like `tune::tune_bayes()`. Starts from `initial` space-filling designs,
then proposes one design per iteration by expected improvement under a
Gaussian process fitted to design means (with Monte Carlo standard
errors as noise). With `constraints`, improvement is weighted by the
modeled probability that each constrained metric meets its bound, so
"the cheapest design with at least 80% power" is answered directly:
`metrics = sim_metric_set(cost, power)`,
`constraints = list(power = 0.8)`.

## Usage

``` r
tune_bayes(object, ...)

# S3 method for class 'sim_spec'
tune_bayes(
  object,
  resamples = NULL,
  ...,
  metrics = NULL,
  param_info = NULL,
  initial = 6,
  iter = 10,
  metric = NULL,
  constraints = NULL,
  no_improve = Inf,
  candidates = 2000,
  control = control_sim(),
  seed = 1
)
```

## Arguments

- object:

  A simulation spec with
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html)
  arguments.

- ...:

  Not used.

- resamples:

  Replicates from
  [`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md),
  or a number of replicates.

- metrics:

  A
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md);
  defaults to the spec's own.

- param_info:

  A dials `parameters` object with changed ranges, from
  `extract_parameter_set_dials(spec) |> update(...)`.

- initial:

  Number of starting designs, or a data frame of them.

- iter:

  Number of proposals after the initial designs.

- metric:

  The metric to race on; defaults to the first in the set.

- constraints:

  Named list of bounds on other metrics' means, e.g.
  `list(power = 0.8)`. A bare number is a floor for metrics to maximize
  and a ceiling for metrics to minimize; `list(cost = c("<=", 5000))` is
  explicit.

- no_improve:

  Stop after this many iterations without a new best design.

- candidates:

  Size of the random candidate pool scored at each iteration.

- control:

  A
  [`control_sim()`](https://jhgille2.github.io/tidyagsim/reference/control_sim.md)
  object.

- seed:

  Seed for generating the space-filling grid.

## Value

A `sim_results` object; `x$trace` lists each proposal.
