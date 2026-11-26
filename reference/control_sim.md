# Control the running of a simulation study

Like `tune::control_grid()`.

## Usage

``` r
control_sim(
  verbose = FALSE,
  parallel = FALSE,
  cache = NULL,
  stop_on_error = FALSE,
  save_outcomes = FALSE,
  shard = NULL
)
```

## Arguments

- verbose:

  Print progress.

- parallel:

  Run design batches in parallel with the future framework
  (`future.apply` must be installed and a plan set, e.g.
  `future::plan(future::multisession, workers = 4)`; this works on
  Windows).

- cache:

  A directory where each finished replicate is saved. Rerunning the same
  study skips them, so an interrupted run resumes where it stopped.

- stop_on_error:

  Stop at the first failing simulation instead of recording a note and
  carrying on.

- save_outcomes:

  Keep every outcome the simulator reports (not only the metrics) in a
  list column `.outcome`; see
  [`collect_outcomes()`](https://jhgille2.github.io/tidyagsim/reference/collect_outcomes.md).

- shard:

  `NULL`, or `c(task, n_tasks)` to run only every `n_tasks`-th
  design/scenario batch, starting at `task`. Used with `cache` to split
  a
  [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md)
  across a Slurm job array;
  [`sim_shard()`](https://jhgille2.github.io/tidyagsim/reference/sim_shard.md)
  reads the task from the environment. A run without `shard` then
  collects everything from the cache.

## Value

A `control_sim` object.

## Examples

``` r
control_sim(verbose = TRUE, cache = tempfile("study-cache"))
#> $verbose
#> [1] TRUE
#> 
#> $parallel
#> [1] FALSE
#> 
#> $cache
#> [1] "/tmp/RtmpSuHJhq/study-cache18d03d8fdf1d"
#> 
#> $stop_on_error
#> [1] FALSE
#> 
#> $save_outcomes
#> [1] FALSE
#> 
#> $shard
#> NULL
#> 
#> attr(,"class")
#> [1] "control_sim"
```
