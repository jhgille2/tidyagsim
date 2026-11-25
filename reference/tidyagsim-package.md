# tidyagsim: tidymodels-style tuning and replication for agricultural simulations

tidyagsim maps the tidymodels grammar onto stochastic simulation
studies:

## Details

|  |  |
|----|----|
| tidymodels | tidyagsim |
| parsnip model spec, [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) | a simulation spec such as [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md), [`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md) or [`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md), with [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) arguments |
| dials parameters | [`extract_parameter_set_dials()`](https://hardhat.tidymodels.org/reference/hardhat-extract.html) and parameter functions such as [`n_lines()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md), [`backcross_generations()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md) |
| rsample resamples | paired replicate seeds and assumption scenarios: [`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md), [`sim_scenarios()`](https://jhgille2.github.io/tidyagsim/reference/sim_scenarios.md) |
| yardstick `metric_set()` | [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md) with directions (maximize power, minimize cost) |
| [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md), `finetune::tune_race_anova()`, [`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md) | [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md), [`tune_race()`](https://jhgille2.github.io/tidyagsim/reference/tune_race.md), [`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md) |
| [`collect_metrics()`](https://jhgille2.github.io/tidyagsim/reference/collect_metrics.md), [`show_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md), [`select_by_one_std_err()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md) | the same names, plus [`pareto_front()`](https://jhgille2.github.io/tidyagsim/reference/pareto_front.md) and [`compare_designs()`](https://jhgille2.github.io/tidyagsim/reference/compare_designs.md) |
| workflowsets | [`sim_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_set.md), [`sim_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_map.md), [`rank_results()`](https://jhgille2.github.io/tidyagsim/reference/rank_results.md) |

## See also

Useful links:

- <https://github.com/jhgille2/tidyagsim>

- <https://jhgille2.github.io/tidyagsim/>

- Report bugs at <https://github.com/jhgille2/tidyagsim/issues>

## Author

**Maintainer**: Jay Gillenwater
<55106933+jhgille2@users.noreply.github.com>
