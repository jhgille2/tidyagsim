# tidyagsim

<!-- badges: start -->
[![R-CMD-check](https://github.com/jhgille2/tidyagsim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jhgille2/tidyagsim/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

tidyagsim brings the grammar of [tidymodels](https://www.tidymodels.org/) to
stochastic simulation studies in plant breeding. Mark the arguments of a
simulation with `tune()`, describe replicate seeds and assumption scenarios,
choose a metric set, and let grid search, racing or Bayesian optimization run
the experiment, without enumerating every combination by hand.

```r
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
|---|---|
| parsnip model spec + `tune()` | `mapping_population()`, `breeding_program()`, `sim_function()` + `tune()` |
| dials | `extract_parameter_set_dials()`, `n_lines()`, `filial_generation()`, ... |
| rsample | `sim_replicates()` (paired seeds), `sim_scenarios()` |
| yardstick `metric_set()` | `sim_metric_set()` |
| `tune_grid()`, `finetune::tune_race_anova()`, `tune_bayes()` | `tune_grid()`, `tune_race()`, `tune_bayes()` (with constraints) |
| `collect_metrics()`, `show_best()`, `select_by_one_std_err()` | same names, plus `pareto_front()`, `compare_designs()` |
| workflowsets | `sim_set()`, `sim_map()`, `rank_results()` |

## Simulators

* **`mapping_population()`**: a fast simulator for QTL mapping populations
  (F2, recombinant inbred lines, BC*s*F*t*, nested association mapping,
  advanced-backcross NAM and 4- or 8-way MAGIC), with genome scans,
  permutation thresholds, costs and time to results. See
  `vignette("mapping-design")`.
* **`breeding_program()`**: a breedsim breeding-program flowchart run by the
  AlphaSimR engine, tuned by `"stage.field"` paths.
* **`sim_function()`**: any R function of its arguments and a seed.

## Installation

```r
# install.packages("pak")
pak::pak("jhgille2/tidyagsim")
```

On Windows this needs only R (4.1 or later). For `breeding_program()` you also
need Python with the breedsim package and the R packages AlphaSimR and
jsonlite. From PowerShell, in a clone of plant_breeding_simulator_flowchart:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e .
```

Then in R, point reticulate at that Python before loading tidyagsim:

```r
Sys.setenv(RETICULATE_PYTHON = "C:/path/to/plant_breeding_simulator_flowchart/.venv/Scripts/python.exe")
Sys.setenv(BREEDSIM_RSCRIPT = "C:/Program Files/R/R-4.5.1/bin/Rscript.exe")  # if Rscript is not on PATH
```

## Documentation

The pkgdown site is at <https://jhgille2.github.io/tidyagsim/>.
