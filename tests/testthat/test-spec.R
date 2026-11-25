test_that("tune() arguments become parameters with default ranges", {
  spec <- mapping_population(type = "RIL", lines = tune(), generation = tune("t"))
  expect_true(is_tune(spec$args$lines))
  args <- tune_args(spec)
  expect_equal(args$name, c("lines", "generation"))
  expect_equal(args$id, c("lines", "t"))
  params <- extract_parameter_set_dials(spec)
  expect_s3_class(params, "parameters")
  expect_equal(params$id, c("lines", "t"))
  narrowed <- update(params, lines = n_lines(c(100, 300)))
  expect_equal(unname(unlist(narrowed$object[[1]]$range)), c(100, 300))
})

test_that("set_args, finalize_sim and unknown arguments", {
  spec <- mapping_population(type = "F2") |> set_args(lines = tune())
  expect_true(is_tune(spec$args$lines))
  expect_error(set_args(spec, nope = 1), "no argument")
  final <- finalize_sim(spec, tibble::tibble(lines = 150L))
  expect_equal(final$args$lines, 150L)
  expect_equal(nrow(tune_args(final)), 0)
})

test_that("shared ids link arguments", {
  spec <- mapping_population(locations = tune("field"), reps = tune("field"))
  expect_equal(extract_parameter_set_dials(spec)$id, "field")
  args <- resolve_args(spec, list(field = 2L))
  expect_equal(c(args$locations, args$reps), c(2L, 2L))
})

test_that("sim_function needs ranges for arguments without defaults", {
  spec <- sim_function(function(a, seed) a, a = tune())
  expect_error(extract_parameter_set_dials(spec), "No default range")
  expect_error(sim_function(function(a) a), "seed")
})

test_that("printing works", {
  expect_output(print(mapping_population(lines = tune())), "Tuned arguments")
})
