# The job-array shard of this process

Reads the Slurm array task (`SLURM_ARRAY_TASK_ID`) and the number of
tasks (`TIDYAGSIM_TASKS`, set by
[`write_slurm_array()`](https://jhgille2.github.io/tidyagsim/reference/write_slurm_array.md),
or `SLURM_ARRAY_TASK_COUNT`). Outside an array it returns `NULL`, so the
same script both runs a shard and, run once more without an array,
collects the finished study.

## Usage

``` r
sim_shard()
```

## Value

`c(task, n_tasks)` or `NULL`.
