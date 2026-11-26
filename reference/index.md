# Package index

## Simulation specs

What to simulate. Set any argument to
[`tune()`](https://hardhat.tidymodels.org/reference/tune.html).

- [`alphasim_program()`](https://jhgille2.github.io/tidyagsim/reference/alphasim_program.md)
  : A simulation spec from an AlphaSimR script
- [`line_pipeline()`](https://jhgille2.github.io/tidyagsim/reference/line_pipeline.md)
  : A soybean-style inbred line breeding pipeline
- [`mapping_population()`](https://jhgille2.github.io/tidyagsim/reference/mapping_population.md)
  : Mapping population design
- [`sim_function()`](https://jhgille2.github.io/tidyagsim/reference/sim_function.md)
  : A simulation spec from any R function
- [`breeding_program()`](https://jhgille2.github.io/tidyagsim/reference/breeding_program.md)
  : A breedsim breeding program as a simulation spec
- [`set_args()`](https://jhgille2.github.io/tidyagsim/reference/set_args.md)
  : Change the arguments of a spec
- [`set_engine()`](https://jhgille2.github.io/tidyagsim/reference/set_engine.md)
  : Change the engine of a spec
- [`finalize_sim()`](https://jhgille2.github.io/tidyagsim/reference/finalize_sim.md)
  : Finalize a spec with chosen values
- [`is_tune()`](https://jhgille2.github.io/tidyagsim/reference/is_tune.md)
  : Is an argument a tune() placeholder?

## AlphaSimR helpers

Founders, and turning a run into breeding metrics.

- [`asr_founders()`](https://jhgille2.github.io/tidyagsim/reference/asr_founders.md)
  : Founder population recipe for AlphaSimR
- [`asr_tracker()`](https://jhgille2.github.io/tidyagsim/reference/asr_tracker.md)
  : Record a breeding program as it runs
- [`asr_summary()`](https://jhgille2.github.io/tidyagsim/reference/asr_summary.md)
  : Summarize a tracked program into breeding metrics
- [`collect_trajectories()`](https://jhgille2.github.io/tidyagsim/reference/collect_trajectories.md)
  : Collect yearly trajectories from a study

## Parameters

dials parameter objects and ranges.

- [`extract_parameter_set_dials(`*`<sim_spec>`*`)`](https://jhgille2.github.io/tidyagsim/reference/extract_parameter_set_dials.sim_spec.md)
  : Extract the parameter set of a spec
- [`n_lines()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`n_families()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`backcross_generations()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`filial_generation()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`population_type()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`magic_founders()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`n_locations()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`n_reps()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`qtl_pve()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`plot_heritability()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  [`generations_per_year()`](https://jhgille2.github.io/tidyagsim/reference/sim_params.md)
  : Parameter objects for simulation studies
- [`sim_param()`](https://jhgille2.github.io/tidyagsim/reference/sim_param.md)
  [`sim_param_values()`](https://jhgille2.github.io/tidyagsim/reference/sim_param.md)
  : A quantitative parameter for any argument

## Replicates and scenarios

- [`sim_replicates()`](https://jhgille2.github.io/tidyagsim/reference/sim_replicates.md)
  : Paired replicate seeds
- [`sim_scenarios()`](https://jhgille2.github.io/tidyagsim/reference/sim_scenarios.md)
  : Assumption scenarios

## Metrics

- [`sim_metric_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric_set.md)
  : Combine metrics into a metric set
- [`sim_metric()`](https://jhgille2.github.io/tidyagsim/reference/sim_metric.md)
  : Define a simulation metric

## Tuning

- [`tune_grid()`](https://jhgille2.github.io/tidyagsim/reference/tune_grid.md)
  : Run every design in a grid
- [`tune_race()`](https://jhgille2.github.io/tidyagsim/reference/tune_race.md)
  : Race designs and drop the clearly worse ones early
- [`tune_bayes()`](https://jhgille2.github.io/tidyagsim/reference/tune_bayes.md)
  : Bayesian optimization of a simulation study
- [`fit_resamples()`](https://jhgille2.github.io/tidyagsim/reference/fit_resamples.md)
  : Run one fixed design on every replicate
- [`control_sim()`](https://jhgille2.github.io/tidyagsim/reference/control_sim.md)
  : Control the running of a simulation study

## Running on a cluster

- [`write_slurm_array()`](https://jhgille2.github.io/tidyagsim/reference/write_slurm_array.md)
  : Write Slurm scripts that split a study across a job array
- [`sim_shard()`](https://jhgille2.github.io/tidyagsim/reference/sim_shard.md)
  : The job-array shard of this process

## Results

- [`collect_metrics()`](https://jhgille2.github.io/tidyagsim/reference/collect_metrics.md)
  : Collect tidy metrics
- [`show_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md)
  [`select_best()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md)
  [`select_by_one_std_err()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md)
  [`select_by_pct_loss()`](https://jhgille2.github.io/tidyagsim/reference/show_best.md)
  : Find the best designs
- [`pareto_front()`](https://jhgille2.github.io/tidyagsim/reference/pareto_front.md)
  : Pareto front of designs
- [`compare_designs()`](https://jhgille2.github.io/tidyagsim/reference/compare_designs.md)
  : Paired comparison of two designs
- [`collect_notes()`](https://jhgille2.github.io/tidyagsim/reference/collect_notes.md)
  : Notes from failed simulations
- [`collect_outcomes()`](https://jhgille2.github.io/tidyagsim/reference/collect_outcomes.md)
  : Collect every outcome the simulator reported
- [`autoplot(`*`<sim_results>`*`)`](https://jhgille2.github.io/tidyagsim/reference/autoplot.sim_results.md)
  [`autoplot(`*`<sim_set>`*`)`](https://jhgille2.github.io/tidyagsim/reference/autoplot.sim_results.md)
  : Plot simulation results

## Comparing several specs

- [`sim_set()`](https://jhgille2.github.io/tidyagsim/reference/sim_set.md)
  : A set of simulation specs to compare
- [`sim_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_map.md)
  : Run a tuner over every spec in a set
- [`rank_results()`](https://jhgille2.github.io/tidyagsim/reference/rank_results.md)
  : Rank designs across a set of specs
- [`extract_sim_result()`](https://jhgille2.github.io/tidyagsim/reference/extract_sim_result.md)
  : Pull one spec's results out of a set

## Backcrossing and pyramiding

Marker-assisted backcrossing, single and multiple introgressions.

- [`backcross_program()`](https://jhgille2.github.io/tidyagsim/reference/backcross_program.md)
  : Marker-assisted backcrossing and gene pyramiding
- [`bc_expected_recovery()`](https://jhgille2.github.io/tidyagsim/reference/bc_expected_recovery.md)
  [`min_population_size()`](https://jhgille2.github.io/tidyagsim/reference/bc_expected_recovery.md)
  : Backcrossing by the numbers

## Mapping population engine

- [`sim_genetic_map()`](https://jhgille2.github.io/tidyagsim/reference/sim_genetic_map.md)
  : Simulate a genetic map
- [`qtl_scan()`](https://jhgille2.github.io/tidyagsim/reference/qtl_scan.md)
  : Genome scan for a single QTL
- [`simulate_mapping_replicate()`](https://jhgille2.github.io/tidyagsim/reference/simulate_mapping_replicate.md)
  : Simulate one mapping population replicate

## Extending tidyagsim

Write a new simulator by implementing these generics.

- [`new_sim_spec()`](https://jhgille2.github.io/tidyagsim/reference/sim_spec.md)
  : Simulation specifications
- [`sim_simulate()`](https://jhgille2.github.io/tidyagsim/reference/sim_simulate.md)
  : Run one design of a spec
- [`sim_check()`](https://jhgille2.github.io/tidyagsim/reference/sim_check.md)
  : Check whether a design is feasible
- [`sim_parameters()`](https://jhgille2.github.io/tidyagsim/reference/sim_parameters.md)
  : Default parameter objects of a spec
- [`tune_args(`*`<sim_spec>`*`)`](https://jhgille2.github.io/tidyagsim/reference/tune_args.sim_spec.md)
  : Tuned arguments of a spec
- [`tunable(`*`<sim_spec>`*`)`](https://jhgille2.github.io/tidyagsim/reference/tunable.sim_spec.md)
  : Parameters that can be tuned, like tidymodels' tunable()
- [`resolve_args()`](https://jhgille2.github.io/tidyagsim/reference/resolve_args.md)
  : Concrete arguments for one design
- [`sim_arg_names()`](https://jhgille2.github.io/tidyagsim/reference/sim_arg_names.md)
  : Argument names a spec accepts
- [`sim_default_metrics()`](https://jhgille2.github.io/tidyagsim/reference/sim_default_metrics.md)
  : Default metric set of a spec
- [`sim_grid()`](https://jhgille2.github.io/tidyagsim/reference/sim_grid.md)
  : Build the candidate grid for a study

## Re-exported

- [`reexports`](https://jhgille2.github.io/tidyagsim/reference/reexports.md)
  [`tune`](https://jhgille2.github.io/tidyagsim/reference/reexports.md)
  [`extract_parameter_set_dials`](https://jhgille2.github.io/tidyagsim/reference/reexports.md)
  [`tunable`](https://jhgille2.github.io/tidyagsim/reference/reexports.md)
  [`tune_args`](https://jhgille2.github.io/tidyagsim/reference/reexports.md)
  [`update`](https://jhgille2.github.io/tidyagsim/reference/reexports.md)
  : Objects exported from other packages
