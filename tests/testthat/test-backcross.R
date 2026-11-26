bc_runs <- function(spec, seeds = 1:30, ...) {
  a <- spec$args
  a[names(list(...))] <- list(...)
  do.call(rbind, lapply(sim_simulate(spec, a, seeds), as.data.frame))
}

test_that("theory helpers", {
  expect_equal(bc_expected_recovery(0:2), c(0.5, 0.75, 0.875))
  expect_equal(min_population_size(1 / 4, 0.99), 17)
  expect_equal(min_population_size(1 / 64, 0.99), 293)
})

test_that("recovery without selection follows 1 - 0.5^(t+1), less the drag", {
  spec <- backcross_program()
  for (bc in c(1, 3)) {
    r <- bc_runs(spec, bc_generations = bc, bg_spacing = Inf, flank_cM = Inf)
    expect_true(all(r$success == 1))
    expect_lt(abs(mean(r$rp_recovery) / 100 - bc_expected_recovery(bc)), 0.04)
    expect_lt(mean(r$rp_recovery) / 100, bc_expected_recovery(bc))   # the carrier chromosome lags
  }
})

test_that("flanking and background selection do their jobs", {
  spec <- backcross_program()
  none <- bc_runs(spec, bg_spacing = Inf, flank_cM = Inf)
  flank <- bc_runs(spec, bg_spacing = Inf, flank_cM = 5)
  both <- bc_runs(spec, flank_cM = 5)
  expect_lt(mean(flank$linkage_drag), mean(none$linkage_drag) / 3)
  expect_gt(mean(both$rp_recovery), mean(flank$rp_recovery) + 3)
  expect_equal(none$datapoints[1], 3 * 100 + 50)                   # foreground only
  expect_gt(mean(both$datapoints), mean(flank$datapoints))
  expect_equal(both$cost, both$plants * 2 + both$datapoints * 1)
})

test_that("more targets need more plants to fix them", {
  spec <- backcross_program(targets = data.frame(chr = c(1, 3, 5), pos = c(50, 30, 90)))
  few <- bc_runs(spec, fix_plants = 20)
  many <- bc_runs(spec, fix_plants = 300)
  expect_lt(mean(few$success), 0.5)
  expect_gt(mean(many$success), 0.85)
  par <- bc_runs(spec, strategy = "parallel", fix_plants = 300)
  expect_equal(par$generations[1], 1 + 3 + 1 + 2 + 1)
  expect_gt(mean(par$rp_recovery), mean(many$rp_recovery))
  expect_gt(mean(par$plants), mean(many$plants))
})

test_that("backcross_program checks its inputs and tunes", {
  expect_error(backcross_program(targets = data.frame(x = 1)), "chr")
  expect_error(backcross_program(n_plantz = 3))
  spec <- backcross_program(targets = data.frame(chr = 30, pos = 1))
  expect_match(sim_check(spec, spec$args), "outside")
  spec <- backcross_program(n_plants = tune(), flank_cM = 5)
  res <- tune_grid(spec, sim_replicates(4, seed = 2), grid = data.frame(n_plants = c(30L, 200L)))
  m <- collect_metrics(res, type = "wide")
  expect_equal(nrow(m), 2)
  expect_true(all(c("rp_recovery", "linkage_drag", "datapoints", "cost") %in% names(m)))
  again <- tune_grid(spec, sim_replicates(4, seed = 2), grid = data.frame(n_plants = c(30L, 200L)))
  expect_equal(collect_metrics(again), collect_metrics(res))
})
