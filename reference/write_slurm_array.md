# Write Slurm scripts that split a study across a job array

Writes `run.sbatch`, an array job in which each task runs `script` with
its shard of the designs, and `collect.sbatch`, which runs after every
task succeeds and assembles the results from the shared cache. The study
script should build its control with
`control_sim(cache = "<shared dir>", shard = sim_shard())` and save the
results (e.g. with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html))
when it is not a shard (`is.null(sim_shard())`).

## Usage

``` r
write_slurm_array(
  script,
  dir,
  tasks,
  cluster = c("ceres", "atlas", "other"),
  account = NULL,
  time = "04:00:00",
  memory = "4G",
  cpus = 1,
  partition = NULL,
  modules = "r",
  container = NULL,
  max_concurrent = NULL
)
```

## Arguments

- script:

  Path of the R script to run (copied into `dir`).

- dir:

  Output directory.

- tasks:

  Number of array tasks.

- cluster:

  `"ceres"`, `"atlas"` or `"other"`.

- account:

  Slurm account (required on Atlas).

- time:

  Walltime per task, e.g. `"04:00:00"`.

- memory:

  Memory per task, e.g. `"8G"`.

- cpus:

  CPUs per task (use with `control_sim(parallel = TRUE)`).

- partition:

  Override the profile's partition.

- modules:

  Environment modules to load, e.g. `"r"`; empty when using a container.

- container:

  Optional Apptainer image with R and tidyagsim; tasks then run
  `apptainer exec <image> Rscript`.

- max_concurrent:

  Limit on simultaneously running tasks.

## Value

The paths written, invisibly.

## Details

Profiles exist for SCINet Ceres and Atlas; other clusters work by
passing `partition` and `account` directly. Submit with
`sbatch run.sbatch`; it submits the collection step itself.

## Examples

``` r
script <- tempfile(fileext = ".R")
writeLines("library(tidyagsim)", script)
write_slurm_array(script, tempfile("slurm"), tasks = 20, cluster = "ceres")
```
