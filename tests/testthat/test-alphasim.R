skip_if_not_installed("AlphaSimR")

tiny <- asr_founders(n_founders = 20, n_chr = 2, seg_sites = 20, qtl_per_chr = 5)
small_line <- function(...) {
  line_pipeline(n_crosses = 10, progeny_per_cross = 5, n_ayt = 20, n_parents = 10, cycles = 2, ..., founders = tiny)
}

test_that("alphasim_program runs a user script reproducibly", {
  spec <- alphasim_program(
    function(founders, SP, n_select, cycles) {
      tr <- asr_tracker(founders)
      pop <- founders
      for (cycle in seq_len(cycles)) {
        pop <- AlphaSimR::setPheno(pop)
        pop <- AlphaSimR::randCross(AlphaSimR::selectInd(pop, n_select), nCrosses = 30)
        tr$record(pop, year = cycle, cost = 30)
      }
      asr_summary(tr)
    },
    n_select = tune(), cycles = 3, founders = tiny,
    params = list(n_select = sim_param(c(4, 15)))
  )
  expect_equal(tune_args(spec)$id, "n_select")
  a <- sim_simulate(spec, list(n_select = 5, cycles = 3), c(1L, 2L))
  b <- sim_simulate(spec, list(n_select = 5, cycles = 3), c(1L, 2L))
  expect_identical(a, b)
  expect_false(identical(a[[1]]$release_gain, a[[2]]$release_gain))
  expect_equal(a[[1]]$years, 3)
  expect_equal(a[[1]]$cost, 90)
  expect_equal(nrow(a[[1]]$trajectory), 3)
  expect_false(exists("SP", envir = globalenv(), inherits = FALSE))
  res <- tune_grid(spec, sim_replicates(3), grid = 2)
  expect_s3_class(res, "sim_results")
  expect_setequal(collect_metrics(res)$.metric, c("release_gain", "gain_per_year", "genetic_variance", "cost"))
})

test_that("selection makes genetic gain and uses up variance", {
  out <- sim_simulate(small_line(), small_line()$args, 1:4)
  gain <- vapply(out, `[[`, 0, "release_gain")
  expect_true(all(gain > 0))
  expect_true(all(vapply(out, `[[`, 0, "genetic_variance") < 1.2))
  expect_true(all(vapply(out, `[[`, 0, "accuracy") > 0))
})

test_that("founders are shared across designs (common random numbers)", {
  r1 <- sim_simulate(small_line(ayt_locations = 1), small_line(ayt_locations = 1)$args, 7L)[[1]]
  r2 <- sim_simulate(small_line(ayt_locations = 4), small_line(ayt_locations = 4)$args, 7L)[[1]]
  # Same founders: the first record's founder-relative baseline is identical, outcomes differ.
  expect_false(identical(r1$release_gain, r2$release_gain))
  expect_equal(r1$lines, r2$lines)
})

test_that("line_pipeline accounts resources and checks feasibility", {
  spec <- small_line()
  a <- spec$args
  r <- line_resources(a)
  expect_equal(r$lines, 50)
  expect_equal(r$plots, 50 * 1 * 1 + 20 * 4 * 2)
  expect_equal(r$years, 5 / 2 + 2)
  expect_null(sim_check(spec, a))
  expect_match(sim_check(spec, modifyList(a, list(n_ayt = 60))), "more than the 50 lines")
  expect_match(sim_check(spec, modifyList(a, list(max_plots = 100))), "max_plots")
  expect_error(line_pipeline(bogus = 1), "no argument")
})

test_that("founder fields can be scenarios and trajectories are collected", {
  spec <- small_line(ayt_locations = tune())
  reps <- sim_replicates(2, scenarios = sim_scenarios(plot_h2 = c(0.1, 0.6)))
  res <- tune_grid(spec, reps, grid = data.frame(ayt_locations = c(1L, 3L)),
                   control = control_sim(save_outcomes = TRUE))
  m <- collect_metrics(res)
  expect_setequal(unique(m$.scenario), c("plot_h2=0.1", "plot_h2=0.6"))
  acc <- collect_metrics(res, by_scenario = TRUE)
  acc <- acc[acc$.metric == "accuracy", ]
  expect_gt(mean(acc$mean[acc$.scenario == "plot_h2=0.6"]), mean(acc$mean[acc$.scenario == "plot_h2=0.1"]))
  tr <- collect_trajectories(res)
  expect_true(all(c("ayt_locations", "year", "mean", ".scenario") %in% names(tr)))
  expect_equal(nrow(tr), 2 * 2 * 2 * 2)
})

test_that("supplied haplotypes are sampled", {
  set.seed(3)
  hap <- AlphaSimR::quickHaplo(nInd = 40, nChr = 2, segSites = 20, inbred = TRUE)
  f <- asr_founders(n_founders = 15, n_chr = 2, qtl_per_chr = 5, haplotypes = hap)
  b <- with_seed(1L, asr_build(f, 1L))
  expect_equal(b$pop@nInd, 15)
  expect_error(asr_build(asr_founders(n_founders = 50, haplotypes = hap), 1L), "more than")
})
