skip_if_no_breedsim <- function() {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("breedsim"), "Python breedsim not installed")
  skip_if(Sys.getenv("BREEDSIM_REPO") == "", "Set BREEDSIM_REPO to the breedsim repository")
}

test_that("breeding programs tune flowchart fields (engine mocked)", {
  skip_if_no_breedsim()
  reticulate::py_run_string("
from unittest.mock import patch
import random
def _fixture(program, panel, settings, rscript):
    reps = next(n.reps for n in program.nodes if n.id == 'pyt')
    return [dict(seed=s, kind='release', year=y, entries=1,
                 gain=-(reps - 2) ** 2 + 0.1 * y + random.Random(s).gauss(0, 0.1))
            for s in panel for y in range(1, program.years + 1)]
_tidyagsim_patch = patch('breedsim.simulation.simulate', side_effect=_fixture)
_tidyagsim_patch.start()
")
  on.exit(reticulate::py_run_string("_tidyagsim_patch.stop()"))
  path <- file.path(Sys.getenv("BREEDSIM_REPO"), "examples", "conventional.json")
  spec <- breeding_program(path, `pyt.reps` = tune(), start_year = 5)
  expect_output(print(spec), "pyt.reps")
  params <- extract_parameter_set_dials(spec)
  expect_equal(params$id, "pyt.reps")
  res <- tune_grid(spec, sim_replicates(3), grid = data.frame(pyt.reps = 1:3),
                   metrics = sim_metric_set(release_gain, gain_per_year, plots))
  expect_equal(select_best(res)$pyt.reps, 2)
  expect_equal(collect_metrics(res, type = "wide")$gain_per_year, rep(0.1, 3), tolerance = 1e-8)
  expect_error(set_args(spec, `nope.reps` = 2), "no argument")
})
