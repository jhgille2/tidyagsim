test_that("tune_grid runs every design on every replicate", {
  res <- tune_grid(toy_spec(), sim_replicates(4), grid = data.frame(x = c(1, 3, 5)), metrics = toy_metrics())
  expect_s3_class(res, "sim_results")
  expect_equal(nrow(res$raw), 12)
  m <- collect_metrics(res)
  expect_equal(nrow(m), 6)
  expect_true(all(c("x", ".metric", "mean", "n", "std_err", ".config") %in% names(m)))
  expect_equal(show_best(res, n = 1)$x, 3)
  expect_equal(select_best(res, metric = "cost")$x, 1)
  wide <- collect_metrics(res, type = "wide")
  expect_true(all(c("y", "y_se", "cost", "cost_se") %in% names(wide)))
  expect_equal(nrow(collect_metrics(res, summarize = FALSE)), 12)
  expect_output(print(res), "Best by y")
})

test_that("generated grids, qualitative parameters and paired seeds", {
  spec <- toy_spec(k = tune())
  res <- tune_grid(spec, sim_replicates(3), grid = 8, metrics = toy_metrics())
  expect_true(nrow(res$configs) >= 6)
  expect_true(all(res$configs$k %in% c("a", "b")))
  # Same seed for the same design reproduces the outcome exactly.
  again <- tune_grid(spec, sim_replicates(3), grid = res$configs[, c("x", "k")], metrics = toy_metrics())
  expect_equal(again$raw$y, res$raw$y)
})

test_that("selection helpers follow tune's rules", {
  res <- tune_grid(toy_spec(noise = 2), sim_replicates(5), grid = data.frame(x = c(2, 2.5, 3, 4)),
                   metrics = toy_metrics())
  simple <- select_by_one_std_err(res, x, metric = "y")
  best <- select_best(res, metric = "y")
  expect_lte(simple$x, best$x)
  expect_true(select_by_pct_loss(res, cost, metric = "y", limit = 100)$x <= best$x)
  front <- pareto_front(res, c("y", "cost"))
  expect_true(nrow(front) >= 1)
  cmp <- compare_designs(res, "Config001", "Config003", metric = "y")
  expect_equal(cmp$n, 5)
  expect_error(show_best(res, metric = "nope"), "not in the metric set")
  constrained <- show_best(res, metric = "y", constraints = list(cost = 25))
  expect_true(all(constrained$x <= 2.5))
})

test_that("scenarios are summarized separately or pooled", {
  spec <- toy_spec()
  r <- sim_replicates(3, scenarios = sim_scenarios(noise = c(0, 1)))
  res <- tune_grid(spec, r, grid = data.frame(x = c(2, 3)), metrics = toy_metrics())
  expect_equal(nrow(res$raw), 12)
  m <- collect_metrics(res)
  expect_true(".scenario" %in% names(m))
  pooled <- collect_metrics(res, by_scenario = FALSE)
  expect_false(".scenario" %in% names(pooled))
  expect_equal(nrow(show_best(res, scenario = "noise=0")), 2)
})

test_that("infeasible designs and failures are recorded, not fatal", {
  spec <- sim_function(function(x, seed) if (x > 4) stop("boom") else list(y = x), x = tune(),
                       params = list(x = sim_param(c(0, 6))),
                       check = function(x) if (x < 1) "too small" else NULL)
  expect_warning(res <- tune_grid(spec, 2, grid = data.frame(x = c(0, 2, 5)),
                                  metrics = sim_metric_set(sim_metric("y"))), "failed")
  expect_equal(res$infeasible$.config, "Config001")
  expect_equal(nrow(res$notes), 1)
  expect_equal(unique(res$raw$.config), "Config002")
})

test_that("cache resumes finished replicates", {
  dir <- tempfile("cache")
  calls <- 0
  counter <- function(x, seed) {
    calls <<- calls + 1
    list(y = x)
  }
  spec <- sim_function(counter, x = tune(), params = list(x = sim_param(c(0, 6))))
  ctrl <- control_sim(cache = dir)
  tune_grid(spec, 3, grid = data.frame(x = 1:2), metrics = sim_metric_set(sim_metric("y")), control = ctrl)
  expect_equal(calls, 6)
  tune_grid(spec, 4, grid = data.frame(x = 1:2), metrics = sim_metric_set(sim_metric("y")), control = ctrl)
  expect_equal(calls, 8)
})

test_that("fit_resamples runs one fixed design", {
  res <- fit_resamples(finalize_sim(toy_spec(), list(x = 3)), 4, metrics = toy_metrics())
  expect_equal(nrow(res$raw), 4)
  expect_error(fit_resamples(toy_spec(), 2, metrics = toy_metrics()), "Fix or finalize")
})

test_that("racing drops clearly worse designs and finishes survivors", {
  res <- tune_race(toy_spec(noise = 0.05), sim_replicates(8), grid = data.frame(x = c(0, 1, 3, 3.1)),
                   metrics = toy_metrics(), burn_in = 3)
  expect_true(nrow(res$trace) >= 1)
  counts <- table(res$raw$.config)
  expect_equal(unname(counts[["Config003"]]), 8)
  expect_lt(unname(counts[["Config001"]]), 8)
  expect_equal(select_best(res)$x, 3)
})

test_that("Bayesian optimization finds the optimum and honours constraints", {
  res <- tune_bayes(toy_spec(noise = 0.05), sim_replicates(2), metrics = toy_metrics(), initial = 4, iter = 6)
  expect_equal(nrow(res$configs), 10)
  expect_lt(abs(select_best(res)$x - 3), 0.75)
  con <- tune_bayes(toy_spec(noise = 0.05), sim_replicates(2), metrics = sim_metric_set(cost, sim_metric("y")),
                    constraints = list(y = -1), initial = 4, iter = 6)
  best <- select_best(con, metric = "cost", constraints = list(y = -1))
  expect_true(best$x >= 1.5 && best$x <= 3.5)
})

test_that("sim sets rank designs across specs", {
  set <- sim_set(low = toy_spec(), high = toy_spec(k = "b"))
  expect_output(print(set), "not run")
  set <- sim_map(set, "tune_grid", resamples = sim_replicates(2), grid = data.frame(x = c(2, 3)),
                 metrics = toy_metrics())
  ranks <- rank_results(set)
  expect_equal(ranks$wflow_id[1], "high")
  expect_equal(nrow(rank_results(set, select_best = TRUE)), 2)
  expect_equal(nrow(collect_metrics(set)), 8)
  expect_s3_class(extract_sim_result(set, "low"), "sim_results")
})

test_that("autoplot draws marginals, fronts and races", {
  skip_if_not_installed("ggplot2")
  res <- tune_grid(toy_spec(k = tune()), sim_replicates(2), grid = 6, metrics = toy_metrics())
  expect_s3_class(ggplot2::autoplot(res), "ggplot")
  expect_s3_class(ggplot2::autoplot(res, type = "pareto"), "ggplot")
  set <- sim_map(sim_set(a = toy_spec(), b = toy_spec(k = "b")), resamples = 2, grid = data.frame(x = 2:3),
                 metrics = toy_metrics())
  expect_s3_class(ggplot2::autoplot(set), "ggplot")
})
