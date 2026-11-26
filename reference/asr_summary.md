# Summarize a tracked program into breeding metrics

Summarize a tracked program into breeding metrics

## Usage

``` r
asr_summary(tracker, ...)
```

## Arguments

- tracker:

  An
  [`asr_tracker()`](https://jhgille2.github.io/tidyagsim/reference/asr_tracker.md).

- ...:

  Extra outcomes to report (for example `plots = 1200`).

## Value

A named list of outcomes: `release_gain` (mean genetic value of the most
recently released lines, or of the last recorded population, in founder
genetic SDs), `gain_per_year`, `genetic_variance` and
`expected_heterozygosity` (last record, relative to the founders),
`cost`, `gain_per_cost` (gain per 1000 currency units), `years`, and
`trajectory`, a data frame of the yearly records.
