#' tidyagsim: tidymodels-style tuning and replication for agricultural simulations
#'
#' tidyagsim maps the tidymodels grammar onto stochastic simulation studies:
#'
#' | tidymodels | tidyagsim |
#' |---|---|
#' | parsnip model spec, `tune()` | a simulation spec such as [mapping_population()], [sim_function()] or [breeding_program()], with `tune()` arguments |
#' | dials parameters | [extract_parameter_set_dials()] and parameter functions such as [n_lines()], [backcross_generations()] |
#' | rsample resamples | paired replicate seeds and assumption scenarios: [sim_replicates()], [sim_scenarios()] |
#' | yardstick `metric_set()` | [sim_metric_set()] with directions (maximize power, minimize cost) |
#' | `tune_grid()`, `finetune::tune_race_anova()`, `tune_bayes()` | [tune_grid()], [tune_race()], [tune_bayes()] |
#' | `collect_metrics()`, `show_best()`, `select_by_one_std_err()` | the same names, plus [pareto_front()] and [compare_designs()] |
#' | workflowsets | [sim_set()], [sim_map()], [rank_results()] |
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data := %||% abort warn inform
#' @importFrom stats sd var runif rnorm rbinom quantile setNames
#' @importFrom utils head
## usethis namespace: end
NULL
