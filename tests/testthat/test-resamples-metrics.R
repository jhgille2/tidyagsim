test_that("replicates are paired seeds", {
  r <- sim_replicates(5, seed = 10)
  expect_equal(r$seed, 10:14)
  expect_equal(r$id[1], "Rep1")
  expect_s3_class(r, "sim_rset")
  expect_error(sim_replicates(0))
  expect_output(print(r), "5 paired seeds")
})

test_that("scenarios cross factors", {
  sc <- sim_scenarios(qtl_pve = c(0.02, 0.08), plot_h2 = 0.4)
  expect_equal(nrow(sc), 2)
  expect_equal(sc$.scenario[1], "qtl_pve=0.02, plot_h2=0.4")
  r <- sim_replicates(3, scenarios = sc)
  expect_equal(nrow(attr(r, "scenarios")), 2)
  expect_equal(nrow(attr(r[1:2, ], "scenarios")), 2)
})

test_that("metric sets take bare names and custom metrics", {
  per_line <- sim_metric("cost_per_line", "minimize", fn = function(o) o$cost / o$lines)
  ms <- sim_metric_set(power, "cost", per_line)
  expect_equal(names(ms), c("power", "cost", "cost_per_line"))
  expect_equal(ms$power$direction, "maximize")
  expect_equal(ms$cost$direction, "minimize")
  expect_equal(tidyagsim:::metric_value(ms$cost_per_line, list(cost = 100, lines = 4)), 25)
  expect_error(sim_metric_set(not_a_metric), "Unknown metric")
  expect_output(print(ms), "power")
})
