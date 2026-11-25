toy_spec <- function(...) {
  sim_function(
    function(x, k = "a", noise = 0.1, seed) {
      list(y = -(x - 3)^2 + (k == "b") + stats::rnorm(1, sd = noise), cost = 10 * x)
    },
    x = tune(), ..., params = list(x = sim_param(c(0, 6), type = "double"), k = sim_param_values(c("a", "b")))
  )
}
toy_metrics <- function() sim_metric_set(sim_metric("y", "maximize"), cost)
small_mapping <- function(...) {
  args <- utils::modifyList(list(n_chromosomes = 3, chromosome_length = 100, marker_spacing = 2, n_perm = 20,
                                 n_background = 100), list(...))
  do.call(mapping_population, args)
}
