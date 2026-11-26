# tidyagsim 0.1.0

* First version: simulation specs with `tune()` arguments, dials parameters,
  paired replicate seeds and assumption scenarios, metric sets, `tune_grid()`,
  `tune_race()`, `tune_bayes()` with constraints, sim sets, tidy results and
  selection helpers, `autoplot()` methods.
* `mapping_population()`: F2, RIL, BCsFt, NAM, AB-NAM and MAGIC populations
  with marker, family and founder-haplotype scans.
* `breeding_program()`: breedsim flowchart programs via reticulate.
* AlphaSimR engine: `alphasim_program()` turns any AlphaSimR script into a
  tunable spec, with reproducible founders shared across designs
  (`asr_founders()`), `asr_tracker()` / `asr_summary()` for breeding metrics
  and `collect_trajectories()` for yearly records. `line_pipeline()` is a
  ready-made inbred line pipeline with G×E, costs and time.
* A data-frame `grid` now supplies ranges for tuned arguments that have no
  default range.
* `line_pipeline()` gains `pyt_plot_cost`, `pyt_error` (a cheaper or noisier
  PYT screen) and `offseason_cost` (winter nursery and greenhouse generations).
* More articles: glossary, cheat sheet, selection intensity versus diversity,
  winter nursery, phenotyping technology break-even, genomic versus phenotypic
  selection.
* Beginner articles (why simulate, how many replicates, choosing a design)
  and case studies (yield trial allocation, gene pyramiding).
