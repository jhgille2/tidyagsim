mapping_defaults <- function() {
  list(
    type = "RIL", lines = 200L, families = 10L, backcrosses = 1L, generation = 6L, founders = "4",
    locations = 1L, reps = 1L,
    qtl_pve = 0.05, plot_h2 = 0.5, gxe_ratio = 0.5, n_background = 1000L, qtl_donor_freq = 0.5,
    n_chromosomes = 20L, chromosome_length = 115, marker_spacing = 1,
    scan_model = "auto", n_perm = 100L, alpha = 0.05, threshold = NULL, window = 10, scan_step = NULL,
    cross_cost = 25, ssd_cost = 2, genotype_cost = 20, plot_cost = 15, generations_per_year = 2L,
    max_plots = Inf, max_cost = Inf, max_years = Inf
  )
}

#' Mapping population design
#'
#' A simulation spec for QTL mapping populations, so their design can be tuned
#' on its own or later as one stage of a breeding program. Each replicate draws
#' inbred founders, places one QTL of stated effect at a random position, builds
#' the population by simulated meiosis, phenotypes the lines, scans the genome
#' and records whether the QTL was detected.
#'
#' @section Population types:
#' Following R/qtl's BC*s*F*t* notation, with *t* the filial generation:
#'
#' | `type` | Construction | Arguments used |
#' |---|---|---|
#' | `"F2"` | Two founders -> F1 -> one selfing (phenotyped as F2:3 family means) | `lines` |
#' | `"RIL"` | Two founders -> F1 -> single-seed descent to F*t* | `lines`, `generation` |
#' | `"BCsFt"` | Two founders -> F1 -> *s* backcrosses to founder 1 -> self to BC*s*F*t* | `lines`, `backcrosses`, `generation` |
#' | `"NAM"` | One common founder x `families` donors -> RILs in each family | `lines`, `families`, `generation` |
#' | `"AB-NAM"` | As NAM, with *s* backcrosses to the common parent in each family | `lines`, `families`, `backcrosses`, `generation` |
#' | `"MAGIC"` | 4 or 8 founders intermated in a funnel -> self to F*t* | `lines`, `founders`, `generation` |
#'
#' `lines` is always the total number of lines, so designs with the same
#' `lines`, `locations` and `reps` use the same field plots and compare fairly.
#' F1 is generation 1: F2 is `generation = 2`, an F6 RIL is `generation = 6`,
#' and BC1F3 is `backcrosses = 1, generation = 3` (one backcross, then two selfings),
#' matching the convention of breedsim's `mapping-program`.
#'
#' @section QTL effect and genetic model:
#' Traits are additive. `qtl_pve` is the share of single-plot phenotypic
#' variance the QTL would explain in a biparental RIL in which it segregates,
#' so its allele effect is the same in every design and realized variance
#' explained differs by design (it is reported as `realized_pve`). In NAM and
#' AB-NAM, donors carry the QTL allele with probability `qtl_donor_freq`; in
#' MAGIC every founder does. A polygenic background of `n_background` loci
#' gives the plot-basis heritability `plot_h2`, and genotype-by-location
#' variance is `gxe_ratio` times the genetic variance, so more locations and
#' reps sharpen line means.
#'
#' @section Detection:
#' A replicate counts as a detection when the LOD within `window` cM of the
#' QTL exceeds the genome-wide threshold (permutations, or `threshold`). A
#' threshold exceedance on any other chromosome is a false positive. See
#' [qtl_scan()] for the models.
#'
#' @section Costs and time:
#' `cost` adds crosses, single-seed-descent plant generations, genotyping and
#' field plots at the unit costs given (placeholders: replace them with your
#' own). `years` is generations divided by `generations_per_year`, rounded up,
#' plus one year for the field trial. `max_plots`, `max_cost` and `max_years`
#' make designs over budget infeasible, so tuners skip them.
#'
#' @param type Population type (see above).
#' @param lines Total lines phenotyped and genotyped.
#' @param families Number of families (NAM and AB-NAM).
#' @param backcrosses Backcross generations *s* (BCsFt and AB-NAM).
#' @param generation Filial generation *t* at phenotyping (RIL, BCsFt, NAM, AB-NAM, MAGIC).
#' @param founders MAGIC founders, `"4"` or `"8"`.
#' @param locations,reps Field locations and replicates per location.
#' @param qtl_pve Variance explained by the QTL in a reference RIL.
#' @param plot_h2 Plot-basis heritability of the polygenic background.
#' @param gxe_ratio Genotype-by-location variance as a multiple of genetic variance.
#' @param n_background Number of background loci.
#' @param qtl_donor_freq Probability a donor or MAGIC founder carries the QTL allele.
#' @param n_chromosomes,chromosome_length,marker_spacing Genome (see [sim_genetic_map()]).
#' @param scan_model `"auto"`, `"marker"`, `"family_marker"` or `"haplotype"`.
#' @param n_perm Permutations per replicate for the threshold.
#' @param alpha Genome-wide type I error rate.
#' @param threshold Fixed LOD threshold (skips permutations; much faster).
#' @param window Distance in cM within which a peak counts as the QTL.
#' @param scan_step Scan every `scan_step` cM instead of every marker.
#' @param cross_cost,ssd_cost,genotype_cost,plot_cost Unit costs.
#' @param generations_per_year Generations per year (field, winter nursery, greenhouse).
#' @param max_plots,max_cost,max_years Budget limits; designs over them are infeasible.
#' @return A `mapping_population` simulation spec.
#' @examples
#' spec <- mapping_population(type = "RIL", lines = tune(), generation = tune(), n_chromosomes = 4)
#' spec
#' extract_parameter_set_dials(spec)
#' @export
mapping_population <- function(type = "RIL", lines = 200L, families = 10L, backcrosses = 1L, generation = 6L,
                               founders = "4", locations = 1L, reps = 1L, qtl_pve = 0.05, plot_h2 = 0.5,
                               gxe_ratio = 0.5, n_background = 1000L, qtl_donor_freq = 0.5, n_chromosomes = 20L,
                               chromosome_length = 115, marker_spacing = 1, scan_model = "auto", n_perm = 100L,
                               alpha = 0.05, threshold = NULL, window = 10, scan_step = NULL, cross_cost = 25,
                               ssd_cost = 2, genotype_cost = 20, plot_cost = 15, generations_per_year = 2L,
                               max_plots = Inf, max_cost = Inf, max_years = Inf) {
  args <- as.list(environment())
  args <- args[names(mapping_defaults())]
  new_sim_spec("mapping_population", args, engine = "tidyagsim")
}

#' @export
sim_arg_names.mapping_population <- function(object) names(mapping_defaults())

#' @export
sim_parameters.mapping_population <- function(object, ...) {
  list(type = population_type(), lines = n_lines(), families = n_families(),
       backcrosses = backcross_generations(), generation = filial_generation(), founders = magic_founders(),
       locations = n_locations(), reps = n_reps(), qtl_pve = qtl_pve(), plot_h2 = plot_heritability(),
       generations_per_year = generations_per_year(),
       n_perm = sim_param(c(20L, 1000L), label = "Permutations"),
       window = sim_param(c(2, 20), type = "double", label = "Detection window (cM)"))
}

#' @export
sim_default_metrics.mapping_population <- function(object) {
  sim_metric_set("power", "false_positive", "location_error", "lines", "cost", "years")
}

#' @export
sim_check.mapping_population <- function(object, args, ...) {
  a <- args
  types <- c("F2", "RIL", "BCsFt", "NAM", "AB-NAM", "MAGIC")
  if (!a$type %in% types) return(sprintf("Unknown population type %s.", a$type))
  d <- resolve_design(a)
  if (d$lines < 10) return("Use at least 10 lines.")
  if (d$s == 0 && d$t < 2) return("Without backcrosses, the filial generation must be at least 2 (F2).")
  if (d$type == "BCsFt" && d$s < 1) return("BCsFt needs at least one backcross; use RIL or F2 for none.")
  if (d$type == "AB-NAM" && d$s < 1) return("AB-NAM needs at least one backcross; use NAM for none.")
  if (d$type == "MAGIC" && !d$founders %in% c(4L, 8L)) return("MAGIC uses 4 or 8 founders.")
  if (d$type %in% c("NAM", "AB-NAM")) {
    if (d$families < 2) return("NAM needs at least two families.")
    if (d$lines / d$families < 5) return(sprintf("%d lines are too few for %d families (fewer than 5 per family).",
                                                  d$lines, d$families))
  }
  plots <- d$lines * a$locations * a$reps
  if (plots > a$max_plots) return(sprintf("%d plots exceed the limit of %s.", plots, format(a$max_plots)))
  info <- design_summary(d$type, d$lines, d$families, d$s, d$t, d$founders)
  years <- ceiling(info$generations / a$generations_per_year) + 1
  if (years > a$max_years) return(sprintf("%d years exceed the limit of %s.", years, format(a$max_years)))
  cost <- info$crosses * a$cross_cost + info$plant_generations * a$ssd_cost + d$lines * a$genotype_cost +
    plots * a$plot_cost
  if (cost > a$max_cost) return(sprintf("Cost %s exceeds the limit of %s.", format(cost), format(a$max_cost)))
  NULL
}

#' @export
sim_simulate.mapping_population <- function(object, args, seeds, ...) {
  lapply(seeds, function(seed) {
    out <- with_seed(seed, simulate_mapping_replicate(args, seed))
    out
  })
}
