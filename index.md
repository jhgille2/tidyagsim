# tidyagsim

tidyagsim brings the grammar of
[tidymodels](https://www.tidymodels.org/) to stochastic simulation
studies in plant breeding. Mark the arguments of a simulation with
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html), describe
replicate seeds and assumption scenarios, choose a metric set, and let
grid search, racing or Bayesian optimization run the experiment, without
enumerating every combination by hand.

``` r

library(tidyagsim)

spec <- mapping_population(type = "RIL", lines = tune(), locations = tune(), qtl_pve = 0.05)

res <- tune_grid(
  spec,
  resamples = sim_replicates(40, scenarios = sim_scenarios(qtl_pve = c(0.03, 0.08))),
  grid = expand.grid(lines = c(100, 200, 400), locations = 1:3),
  metrics = sim_metric_set(power, location_error, cost)
)

collect_metrics(res)                                  # tidy means and Monte Carlo SEs
select_by_one_std_err(res, cost, metric = "power")    # cheapest design near the best
autoplot(res)
```

| tidymodels | tidyagsim |
|----|----|
| parsnip model spec + [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) | [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md), [`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md), [`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md) + [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) |
| dials | [`extract_parameter_set_dials()`](https://hardhat.tidymodels.org/reference/hardhat-extract.html), [`n_lines()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md), [`filial_generation()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md), … |
| rsample | [`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md) (paired seeds), [`sim_scenarios()`](https://jhgille2.github.io/tidyagsim/reference/sim_scenarios.md) |
| yardstick `metric_set()` | [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md) |
| [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md), `finetune::tune_race_anova()`, [`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md) | [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md), [`tune_race()`](https://jhgille2.github.io/tidyagsim/reference/tune_race.md), [`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md) (with constraints) |
| [`collect_metrics()`](https://jhgille2.github.io/tidyagsim/reference/collect_metrics.md), [`show_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md), [`select_by_one_std_err()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md) | same names, plus [`pareto_front()`](https://jhgille2.github.io/tidyagsim/reference/pareto_front.md), [`compare_designs()`](https://jhgille2.github.io/tidyagsim/reference/compare_designs.md) |
| workflowsets | [`sim_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_set.md), [`sim_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_map.md), [`rank_results()`](https://jhgille2.github.io/tidyagsim/reference/rank_results.md) |

## Simulators

- **[`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md)**:
  a fast simulator for QTL mapping populations (F2, recombinant inbred
  lines, BC*s*F*t*, nested association mapping, advanced-backcross NAM
  and 4- or 8-way MAGIC), with genome scans, permutation thresholds,
  costs and time to results. See
  [`vignette("mapping-design")`](https://jhgille2.github.io/tidyagsim/articles/mapping-design.md).
- **[`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md)**:
  a breedsim breeding-program flowchart run by the AlphaSimR engine,
  tuned by `"stage.field"` paths.
- **[`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md)**:
  any R function of its arguments and a seed.

## Installation

``` r

# install.packages("pak")
pak::pak("jhgille2/tidyagsim")
```

On Windows this needs only R (4.1 or later). For
[`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md)
you also need Python with the breedsim package and the R packages
AlphaSimR and jsonlite. From PowerShell, in a clone of
plant_breeding_simulator_flowchart:

``` powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e .
```

Then in R, point reticulate at that Python before loading tidyagsim:

``` r

Sys.setenv(RETICULATE_PYTHON = "C:/path/to/plant_breeding_simulator_flowchart/.venv/Scripts/python.exe")
Sys.setenv(BREEDSIM_RSCRIPT = "C:/Program Files/R/R-4.5.1/bin/Rscript.exe")  # if Rscript is not on PATH
```

## Documentation

The pkgdown site is at <https://jhgille2.github.io/tidyagsim/>.
