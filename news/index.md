# Changelog

## tidyagsim 0.1.0

- First version: simulation specs with
  [`tune()`](https://hardhat.tidymodels.org/reference/tune.html)
  arguments, dials parameters, paired replicate seeds and assumption
  scenarios, metric sets,
  [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md),
  [`tune_race()`](https://jhgille2.github.io/tidyagsim/reference/tune_race.md),
  [`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md)
  with constraints, sim sets, tidy results and selection helpers,
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods.
- [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md):
  F2, RIL, BCsFt, NAM, AB-NAM and MAGIC populations with marker, family
  and founder-haplotype scans.
- [`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md):
  breedsim flowchart programs via reticulate.
- AlphaSimR engine:
  [`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md)
  turns any AlphaSimR script into a tunable spec, with reproducible
  founders shared across designs
  ([`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md)),
  [`asr_tracker()`](https://jhgille2.github.io/tidyagsim/reference/asr_tracker.md)
  /
  [`asr_summary()`](https://jhgille2.github.io/tidyagsim/reference/asr_summary.md)
  for breeding metrics and
  [`collect_trajectories()`](https://jhgille2.github.io/tidyagsim/reference/collect_trajectories.md)
  for yearly records.
  [`line_pipeline()`](https://jhgille2.github.io/tidyagsim/reference/line_pipeline.md)
  is a ready-made inbred line pipeline with G×E, costs and time.
- A data-frame `grid` now supplies ranges for tuned arguments that have
  no default range.
- [`line_pipeline()`](https://jhgille2.github.io/tidyagsim/reference/line_pipeline.md)
  gains `pyt_plot_cost`, `pyt_error` (a cheaper or noisier PYT screen)
  and `offseason_cost` (winter nursery and greenhouse generations).
- [`backcross_program()`](https://jhgille2.github.io/tidyagsim/reference/backcross_program.md):
  marker-assisted backcrossing and gene pyramiding (single or multiple
  targets, simultaneous or parallel streams), with
  [`bc_expected_recovery()`](https://jhgille2.github.io/tidyagsim/reference/bc_expected_recovery.md)
  and
  [`min_population_size()`](https://jhgille2.github.io/tidyagsim/reference/bc_expected_recovery.md)
  from backcross theory, and an article on planning a backcross.
- More articles: glossary, cheat sheet, selection intensity versus
  diversity, winter nursery, phenotyping technology break-even, genomic
  versus phenotypic selection.
- Beginner articles (why simulate, how many replicates, choosing a
  design) and case studies (yield trial allocation, gene pyramiding).
