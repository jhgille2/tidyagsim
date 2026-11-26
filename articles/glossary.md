# Glossary

Simulation studies borrow words from statistics, machine learning and
plant breeding. This page defines the ones used across tidyagsim, in
plain language, with the function that goes with each.

## The study

- **Simulation spec**:

  A recipe for one simulated experiment: which simulator to run and with
  what settings. Made with
  [`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md),
  [`line_pipeline()`](https://jhgille2.github.io/tidyagsim/reference/line_pipeline.md),
  [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md)
  or
  [`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md).
  A spec does nothing until a tuner runs it.

- **Design**: One combination of the settings you control, such as “60
  lines at 4 locations with 2 reps”. A study compares designs.

- **[`tune()`](https://hardhat.tidymodels.org/reference/tune.html)**:

  Put [`tune()`](https://hardhat.tidymodels.org/reference/tune.html) in
  place of a setting to say “try several values of this”. Every argument
  you mark becomes a column of the grid of designs.

- **Parameter range**:

  The values a tuned setting may take, as a dials object. See them with
  [`extract_parameter_set_dials()`](https://hardhat.tidymodels.org/reference/hardhat-extract.html),
  change them with [`update()`](https://rdrr.io/r/stats/update.html), or
  simply list the designs you want in a data frame passed as `grid`.

- **Grid**:

  The table of designs to run. You can give it yourself
  ([`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) is handy)
  or let `tune_grid(grid = 20)` spread 20 designs over the ranges.

- **Assumption, scenario**:

  A setting you do not control and are unsure of, such as heritability
  or G×E. Rather than choosing it, you run every design under several
  values of it with
  [`sim_scenarios()`](https://jhgille2.github.io/tidyagsim/reference/sim_scenarios.md),
  and see whether the best design changes.

## Randomness

- **Replicate**:

  One run of a design, with its own random founders, crosses and plot
  errors. A study runs each design on many replicates
  ([`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md))
  and reports averages.

- **Seed**: The number that fixes a replicate’s random draws. The same
  seed always gives the same replicate, so results are reproducible.

- **Common random numbers (paired seeds)**:

  Every design runs on the *same* seeds, so they start from the same
  founders. Differences between designs are then measured within each
  seed, which removes shared luck from the comparison.

![Two designs measured on the same five seeds: lines connect each seed's
pair of values. Values vary a lot between seeds but design B is higher
on every seed.](glossary_files/figure-html/unnamed-chunk-2-1.png)

- **Monte Carlo standard error (SE)**:

  How much an average over replicates would move if you reran the study
  with new seeds. It shrinks with the square root of the number of
  replicates, so 4 times the replicates halves it. Reported by
  [`collect_metrics()`](https://jhgille2.github.io/tidyagsim/reference/collect_metrics.md)
  as `std_err`.

## Measuring designs

- **Metric**:

  One number each replicate reports, such as `release_gain` or `cost`.
  Collected in a
  [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md),
  which also records whether bigger or smaller is better (the metric’s
  *direction*). The first metric is the one tuners optimize by default.

- **Trajectory**:

  A replicate’s yearly record (mean genetic value, variance, cost), kept
  with `control_sim(save_outcomes = TRUE)` and gathered with
  [`collect_trajectories()`](https://jhgille2.github.io/tidyagsim/reference/collect_trajectories.md).

## Choosing

- **Constraint**:

  A limit on a metric’s average, such as
  `constraints = list(cost = 30000)`. Designs that break it are left out
  when ranking.

- **One-standard-error rule**:

  Among designs whose average is within one SE of the best, pick the
  simplest or cheapest
  ([`select_by_one_std_err()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md)).
  It avoids paying for a “best” that is only luck.

- **Pareto front**:

  When goals conflict, such as gain now against diversity for later, the
  designs that no other design beats on every goal at once
  ([`pareto_front()`](https://jhgille2.github.io/tidyagsim/reference/pareto_front.md)).
  Every design off the front can be improved for free.

![Scatter of designs by gain and diversity; the designs on the
upper-right edge are joined as the Pareto
front.](glossary_files/figure-html/unnamed-chunk-3-1.png)

- **Racing**:

  Run every design on a few replicates, drop those clearly worse than
  the leader, and spend further replicates on the rest
  ([`tune_race()`](https://jhgille2.github.io/tidyagsim/reference/tune_race.md)).

- **Bayesian optimization**:

  Fit a statistical model to the designs run so far and use it to choose
  the most promising design to run next
  ([`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md)).
  Useful when each design is expensive and there are too many to try.

## Breeding terms used in the examples

- **Founders**:

  The starting population of a simulation, made fresh for each replicate
  by
  [`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md).

- **Cycle**: One round of crossing, line development, testing and
  selection of new parents.

- **Single-seed descent (SSD)**: Advancing each line by selfing one seed
  per plant per generation until it is nearly inbred, often in a
  greenhouse or winter nursery.

- **PYT, AYT**: Preliminary and advanced yield trials. Many lines with
  few plots, then few lines with many plots.

- **Heritability (plot and entry mean)**:

  How much of the variation among plots is genetic. *Plot* heritability
  is for a single plot. The heritability of a line’s mean over locations
  and reps is higher, which is why trials have several plots.

- **G×E (genotype by environment)**: Lines rank differently in different
  environments. Extra locations average it out, and extra reps at the
  same location do not.

- **Selection intensity**: How hard you select: keeping the top 5% is
  more intense than keeping the top 30%. More intensity gives more gain
  now and uses up variance faster.

- **Accuracy**: The correlation between what you select on (trial means,
  predictions) and true genetic value.

- **Genetic variance, expected heterozygosity**: Two measures of the
  diversity left for future selection, both reported relative to the
  founders.

- **Gain per year, gain per cost**: Genetic gain divided by the years or
  money it took. They are the fair way to compare programs that run at
  different speeds or budgets.

## Backcrossing terms

- **Recurrent parent (RP), donor**: In a backcross, the elite line you
  keep crossing back to is the recurrent parent. The donor supplies the
  gene you want to move in.

- **Foreground, recombinant and background selection**: Marker selection
  in three stages: keep plants carrying the target gene, then plants
  whose flanking markers already come from the RP, then plants with the
  most RP alleles across the genome.

- **Linkage drag**: The donor DNA that travels with the target gene
  because it sits close by. Flanking markers near the target cut it
  down.

- **RP genome recovery**:

  The share of the converted line’s genome that comes from the recurrent
  parent. Without background selection it averages 1 − 0.5^(t+1) after
  *t* backcrosses.

- **Pyramiding**: Combining several target genes in one line. They can
  be brought in together in one backcross stream, or separately in
  parallel streams that are crossed together at the end.

- **Introgression, bridging population**: Bringing useful variation from
  unadapted germplasm (such as wild soybean) into an elite program. A
  bridging population crosses exotic and elite lines and selects among
  their progeny before the best go on into the elite crossing block.
