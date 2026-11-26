# tidyagsim

<!-- badges: start -->
[![R-CMD-check](https://github.com/jhgille2/tidyagsim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jhgille2/tidyagsim/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

tidyagsim extends [AlphaSimR](https://cran.r-project.org/package=AlphaSimR)
with the grammar of [tidymodels](https://www.tidymodels.org/). Write a
breeding scheme once, mark the choices you want to compare with `tune()`, and
let tidyagsim run every design on the same replicate seeds, summarize the
results with Monte Carlo standard errors, and help you pick. No nested loops,
no hand-built results tables.

```r
library(tidyagsim)

# A soybean-style line pipeline on AlphaSimR: how should 480 AYT plots be split?
pipeline <- line_pipeline(n_ayt = tune(), ayt_locations = tune(), ayt_reps = tune())

designs <- expand.grid(ayt_locations = c(1, 2, 4, 8), ayt_reps = c(1, 2))
designs$n_ayt <- 480 / (designs$ayt_locations * designs$ayt_reps)

res <- tune_grid(
  pipeline,
  resamples = sim_replicates(20, scenarios = sim_scenarios(gxe_var = c(0.2, 1))),
  grid = designs
)

collect_metrics(res)                          # tidy means and Monte Carlo SEs
show_best(res, metric = "gain_per_year")
autoplot(res)
```

Your own AlphaSimR scripts work the same way:

```r
library(AlphaSimR)

recurrent <- alphasim_program(
  function(founders, SP, n_select, cycles) {
    tracker <- asr_tracker(founders)
    pop <- founders
    for (cycle in seq_len(cycles)) {
      pop <- randCross(selectInd(setPheno(pop), n_select), nCrosses = 200)
      tracker$record(pop, year = cycle, cost = 200 * 15)
    }
    asr_summary(tracker)
  },
  n_select = tune(), cycles = 10,
  params = list(n_select = sim_param(c(4, 60)))
)
tune_bayes(recurrent, sim_replicates(10), constraints = list(genetic_variance = 0.2))
```

| tidymodels | tidyagsim |
|---|---|
| parsnip model spec + `tune()` | `alphasim_program()`, `line_pipeline()`, `mapping_population()`, `sim_function()` + `tune()` |
| dials | `extract_parameter_set_dials()`, `n_lines()`, `filial_generation()`, ... |
| rsample | `sim_replicates()` (paired seeds), `sim_scenarios()` |
| yardstick `metric_set()` | `sim_metric_set()` |
| `tune_grid()`, `finetune::tune_race_anova()`, `tune_bayes()` | `tune_grid()`, `tune_race()`, `tune_bayes()` (with constraints) |
| `collect_metrics()`, `show_best()`, `select_by_one_std_err()` | same names, plus `pareto_front()`, `compare_designs()` |
| workflowsets | `sim_set()`, `sim_map()`, `rank_results()` |

## Simulators

* **`alphasim_program()`**: any AlphaSimR breeding scheme, written as a
  function of the founders, `SP` and your design arguments. tidyagsim builds
  reproducible founders for each replicate (`asr_founders()`), and
  `asr_tracker()` / `asr_summary()` turn a run into gain, gain per year, gain
  per cost, remaining genetic variance and diversity.
* **`line_pipeline()`**: a ready-made inbred line pipeline (crossing, single-seed
  descent, PYT, AYT, release and parent recycling) with G×E, costs and time.
* **`mapping_population()`**: a fast simulator for QTL mapping populations
  (F2, recombinant inbred lines, BC*s*F*t*, nested association mapping,
  advanced-backcross NAM and 4- or 8-way MAGIC), with genome scans,
  permutation thresholds, costs and time to results.
* **`backcross_program()`**: marker-assisted backcrossing of one or more
  target genes, with foreground, flanking-marker and background selection,
  simultaneous or parallel (pyramiding) strategies, and recurrent parent
  recovery, linkage drag, time and genotyping cost as metrics.
* **`sim_function()`**: any R function of its arguments and a seed.
* **`breeding_program()`**: a breedsim flowchart (from the desktop app),
  tuned by `"stage.field"` paths.

## Learning tidyagsim

New to breeding simulation? The articles build up one idea at a time:

1. [Getting started](https://jhgille2.github.io/tidyagsim/articles/tidyagsim.html): the grammar in five minutes.
2. [Why simulate?](https://jhgille2.github.io/tidyagsim/articles/why-simulate.html): checking a simulation against the breeder's equation.
3. [How many replicates?](https://jhgille2.github.io/tidyagsim/articles/replicates.html): Monte Carlo error and paired seeds.
4. [Choosing a design](https://jhgille2.github.io/tidyagsim/articles/choosing-designs.html): budgets, the one-standard-error rule and Pareto fronts.

Then the case studies: [yield trial allocation](https://jhgille2.github.io/tidyagsim/articles/yield-trials.html),
[gene pyramiding](https://jhgille2.github.io/tidyagsim/articles/gene-pyramiding.html)
and [mapping population design](https://jhgille2.github.io/tidyagsim/articles/mapping-design.html).

## Installation

```r
# install.packages("pak")
pak::pak("jhgille2/tidyagsim")
```

On Windows this needs only R (4.1 or later); AlphaSimR installs from CRAN as a
binary. For `breeding_program()` you also need Python with the breedsim
package and the R package jsonlite. From PowerShell, in a clone of plant_breeding_simulator_flowchart:

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
