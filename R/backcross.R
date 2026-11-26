# Marker-assisted backcrossing and pyramiding --------------------------------------------------

backcross_defaults <- function() {
  list(strategy = "simultaneous", bc_generations = 3, n_plants = 100, n_keep = 1, flank_cM = 10,
       bg_spacing = 20, bg_from = 1, fix_plants = 50, pyramid_plants = 50,
       n_chromosomes = 20, chromosome_length = 125, plant_cost = 2, datapoint_cost = 1,
       generations_per_year = 3)
}

#' Marker-assisted backcrossing and gene pyramiding
#'
#' Simulates introgressing one or more target genes from a donor into a
#' recurrent parent (RP) by marker-assisted backcrossing, then fixing them by
#' selfing. Every argument can be tuned, so the classic design questions
#' (how many plants per generation, how many backcross generations, how close
#' the flanking markers should be, how dense the background markers, and
#' whether to introgress several genes together or in parallel streams) become
#' one call to [tune_grid()] or [tune_bayes()].
#'
#' @section Selection in each backcross generation:
#' Each generation grows `n_plants` plants and selects in three stages, as in
#' Hospital and Charcosset (1997) and Frisch et al. (1999):
#' 1. **Foreground:** plants carrying the donor allele at every target, read
#'    from a diagnostic marker on each gene.
#' 2. **Recombinant:** among carriers, those with the most recombinant
#'    flanking markers, placed `flank_cM` on each side of each target. This
#'    limits linkage drag. Use `flank_cM = Inf` to skip this stage.
#' 3. **Background:** from generation `bg_from` on, the remaining plants are
#'    ranked by the share of RP alleles at markers every `bg_spacing` cM. Use
#'    `bg_spacing = Inf` to skip it.
#'
#' The best `n_keep` plants are backcrossed to the RP. After
#' `bc_generations`, the best plant is selfed, `fix_plants` progeny are grown,
#' and a plant homozygous for every target is kept, with ties broken by
#' background.
#'
#' @section Several targets:
#' With more than one target, `strategy = "simultaneous"` backcrosses all of
#' them together, selecting carriers of every target each generation.
#' `strategy = "parallel"` follows Peng et al. (2014a, b). Each target gets its
#' own backcross stream of `n_plants` per generation, and each conversion is
#' fixed by selfing. The conversions are then pyramided by crossing them in
#' pairs, then crossing those hybrids, and so on (a symmetric schedule, as in
#' Ishii and Yonezawa, 2007). Each crossing step grows `pyramid_plants` and
#' keeps the plant carrying every target with the best background. A final
#' selfing fixes the pyramid.
#'
#' @section Outcomes:
#' * `success`: 1 if the final plant is homozygous for the donor allele at
#'   every target, else 0.
#' * `rp_recovery`: percent of the final plant's genome from the RP.
#' * `donor_cM`: donor genome left, in cM per haploid genome.
#' * `linkage_drag`: mean length in cM of the donor segment around each
#'   target.
#' * `generations` and `years`: crossing and selfing generations, and years
#'   at `generations_per_year`. Parallel streams run at the same time.
#' * `plants`, `datapoints` (marker data points) and `cost`.
#'
#' @param ... Arguments to set or tune: `strategy`, `bc_generations`,
#'   `n_plants`, `n_keep`, `flank_cM`, `bg_spacing`, `bg_from`, `fix_plants`,
#'   `pyramid_plants`, `n_chromosomes`, `chromosome_length`, `plant_cost`,
#'   `datapoint_cost`, `generations_per_year`.
#' @param targets Data frame with columns `chr` and `pos` (cM): where the
#'   target genes are.
#' @return A `backcross_program` spec.
#' @references
#' Hospital F, Charcosset A (1997) Marker-assisted introgression of quantitative
#' trait loci. Genetics 147:1469-1485. \doi{10.1093/genetics/147.3.1469}
#'
#' Frisch M, Bohn M, Melchinger AE (1999) Comparison of selection strategies for
#' marker-assisted backcrossing of a gene. Crop Science 39:1295-1301.
#'
#' Hospital F (2001) Size of donor chromosome segments around introgressed loci
#' and reduction of linkage drag in marker-assisted backcross programs.
#' Genetics 158:1363-1379. \doi{10.1093/genetics/158.3.1363}
#'
#' Peng T, Sun X, Mumm RH (2014a) Optimized breeding strategies for multiple
#' trait integration: I. Minimizing linkage drag in single event introgression.
#' Molecular Breeding 33:89-104. \doi{10.1007/s11032-013-9936-7}
#'
#' Peng T, Sun X, Mumm RH (2014b) Optimized breeding strategies for multiple
#' trait integration: II. Process efficiency in event pyramiding and trait
#' fixation. Molecular Breeding 33:105-115. \doi{10.1007/s11032-013-9937-6}
#'
#' Ishii T, Yonezawa K (2007) Optimization of the marker-based procedures for
#' pyramiding genes from multiple donor lines: I. Schedule of crossing between
#' the donor lines. Crop Science 47:537-546.
#' @examples
#' spec <- backcross_program(n_plants = tune(), flank_cM = tune(), n_chromosomes = 5)
#' res <- tune_grid(spec, sim_replicates(5), grid = expand.grid(n_plants = c(50, 150), flank_cM = c(5, 20)))
#' collect_metrics(res, type = "wide")
#' @export
backcross_program <- function(..., targets = data.frame(chr = 1, pos = 50)) {
  args <- list(...)
  targets <- as.data.frame(targets)
  if (!all(c("chr", "pos") %in% names(targets)) || !nrow(targets)) {
    abort("`targets` needs columns `chr` and `pos`.")
  }
  defaults <- backcross_defaults()
  defaults[names(args)] <- args
  spec <- new_sim_spec("backcross_program", defaults, targets = targets[, c("chr", "pos")])
  check_arg_names(spec, names(args))
  spec
}

#' @export
sim_arg_names.backcross_program <- function(object) names(backcross_defaults())

#' @export
sim_parameters.backcross_program <- function(object, ...) {
  list(strategy = sim_param_values(c("simultaneous", "parallel"), label = "Strategy"),
       bc_generations = backcross_generations(c(1L, 6L)),
       n_plants = sim_param(c(20L, 1000L), label = "Plants per generation"),
       n_keep = sim_param(c(1L, 10L), label = "Plants backcrossed"),
       flank_cM = sim_param(c(1, 40), type = "double", label = "Flanking marker distance (cM)"),
       bg_spacing = sim_param(c(5, 60), type = "double", label = "Background marker spacing (cM)"),
       bg_from = sim_param(c(1L, 6L), label = "First generation with background selection"),
       fix_plants = sim_param(c(10L, 500L), label = "Plants grown to fix the targets"),
       pyramid_plants = sim_param(c(10L, 1000L), label = "Plants per pyramiding cross"),
       generations_per_year = generations_per_year(c(1L, 4L)))
}

#' @export
sim_default_metrics.backcross_program <- function(object) {
  sim_metric_set(sim_metric("success", "maximize", label = "Targets fixed"),
                 sim_metric("rp_recovery", "maximize", label = "Recurrent parent genome", units = "%"),
                 sim_metric("linkage_drag", "minimize", label = "Donor segment around targets", units = "cM"),
                 "years", sim_metric("plants", "minimize", label = "Plants grown"),
                 sim_metric("datapoints", "minimize", label = "Marker data points"), "cost")
}

#' @export
sim_check.backcross_program <- function(object, args, ...) {
  a <- args
  if (!a$strategy %in% c("simultaneous", "parallel")) return("strategy is simultaneous or parallel.")
  if (a$n_keep > a$n_plants) return("n_keep is more than n_plants.")
  if (a$bc_generations < 1) return("Use at least one backcross generation.")
  t <- object$targets
  if (any(t$chr > a$n_chromosomes) || any(t$pos < 0 | t$pos > a$chromosome_length)) {
    return("A target lies outside the genome.")
  }
  NULL
}

#' @export
sim_simulate.backcross_program <- function(object, args, seeds, ...) {
  lapply(seeds, function(seed) with_seed(seed, simulate_backcross(args, object$targets)))
}

# The engine ---------------------------------------------------------------------------------

# Genome on a 1 cM grid, with the targets and markers placed on it. Founder 1 is
# the recurrent parent and founder 2 the donor.
bc_genome <- function(a, targets) {
  map <- sim_genetic_map(a$n_chromosomes, a$chromosome_length, 1)
  locus <- function(chr, pos) which(map$chr == chr & map$pos == round(pos))
  tgt <- vapply(seq_len(nrow(targets)), function(i) locus(targets$chr[i], targets$pos[i]), 0L)
  flanks <- if (is.finite(a$flank_cM)) {
    lapply(seq_len(nrow(targets)), function(i) {
      lo <- max(targets$pos[i] - a$flank_cM, 0)
      hi <- min(targets$pos[i] + a$flank_cM, a$chromosome_length)
      c(locus(targets$chr[i], lo), locus(targets$chr[i], hi))
    })
  } else list()
  bg <- if (is.finite(a$bg_spacing)) which(map$pos %% a$bg_spacing == 0) else integer()
  list(map = map, r = recombination_fractions(map$chr, map$pos), targets = tgt, flanks = flanks, bg = bg)
}

bc_counter <- function() {
  env <- new.env()
  env$plants <- 0
  env$datapoints <- 0
  env
}

# Genotype calls: number of donor alleles at loci (rows) for plants (columns).
donor_dose <- function(pop, loci) (pop$h1[loci, , drop = FALSE] == 2L) + (pop$h2[loci, , drop = FALSE] == 2L)

# Pick `n` plants from `pop` that carry every target in `need` (minimum dose), with
# the recombinant and background stages; counts marker data points.
bc_select <- function(pop, g, need, n, counter, use_flanks = TRUE, use_bg = TRUE, min_dose = 1L) {
  k <- ncol(pop$h1)
  counter$datapoints <- counter$datapoints + k * length(need)
  dose <- donor_dose(pop, g$targets[need])
  ok <- which(colSums(dose >= min_dose) == length(need))
  if (!length(ok)) ok <- seq_len(k)                       # target lost: carry on, success will be 0
  cand <- ok
  if (use_flanks && length(g$flanks)) {
    fl <- unlist(g$flanks[need])
    counter$datapoints <- counter$datapoints + length(cand) * length(fl)
    # A flank is "clean" when it carries no donor allele.
    clean <- colSums(donor_dose(list(h1 = pop$h1[, cand, drop = FALSE], h2 = pop$h2[, cand, drop = FALSE]), fl) == 0)
    cand <- cand[clean == max(clean)]
  }
  if (use_bg && length(g$bg) && length(cand) > n) {
    counter$datapoints <- counter$datapoints + length(cand) * length(g$bg)
    sub <- list(h1 = pop$h1[, cand, drop = FALSE], h2 = pop$h2[, cand, drop = FALSE])
    rp <- colSums(donor_dose(sub, g$bg) == 0)
    cand <- cand[order(-rp)]
  }
  # Not enough plants passed: fill with the next best carriers.
  if (length(cand) < n) cand <- c(cand, setdiff(ok, cand))
  if (length(cand) < n) cand <- c(cand, setdiff(seq_len(k), cand))
  cand[seq_len(min(n, length(cand)))]
}

subset_pop <- function(pop, idx) list(h1 = pop$h1[, idx, drop = FALSE], h2 = pop$h2[, idx, drop = FALSE])

# Backcross stream for targets `need`, starting from F1 plants; returns the selected plant(s).
bc_stream <- function(f1, rp, g, need, a, counter) {
  parents <- f1
  for (gen in seq_len(a$bc_generations)) {
    per <- ceiling(a$n_plants / ncol(parents$h1))
    idx <- rep(seq_len(ncol(parents$h1)), each = per)[seq_len(a$n_plants)]
    pop <- list(h1 = gamete(parents, idx, g$r), h2 = rp$h1[, rep(1L, a$n_plants), drop = FALSE])
    counter$plants <- counter$plants + a$n_plants
    keep <- bc_select(pop, g, need, a$n_keep, counter, use_bg = gen >= a$bg_from)
    parents <- subset_pop(pop, keep)
  }
  subset_pop(parents, 1L)
}

# Self one plant and keep a progeny homozygous for every target in `need`.
bc_fix <- function(plant, g, need, n, counter) {
  idx <- rep(1L, n)
  pop <- list(h1 = gamete(plant, idx, g$r), h2 = gamete(plant, idx, g$r))
  counter$plants <- counter$plants + n
  subset_pop(pop, bc_select(pop, g, need, 1L, counter, use_flanks = FALSE, min_dose = 2L))
}

bc_outcomes <- function(plant, g, a, counter, generations) {
  n_t <- length(g$targets)
  dose <- donor_dose(plant, g$targets)
  donor <- (plant$h1 == 2L) + (plant$h2 == 2L)
  # Donor segment around each target on each haplotype that carries it.
  drag <- unlist(lapply(seq_len(n_t), function(i) {
    t <- g$targets[i]
    chr <- g$map$chr == g$map$chr[t]
    vapply(list(plant$h1[, 1], plant$h2[, 1]), function(h) {
      if (h[t] != 2L) return(NA_real_)
      on <- which(chr)
      run <- h[on] == 2L
      pos <- match(t, on)
      lo <- pos
      while (lo > 1 && run[lo - 1]) lo <- lo - 1
      hi <- pos
      while (hi < length(on) && run[hi + 1]) hi <- hi + 1
      g$map$pos[on[hi]] - g$map$pos[on[lo]]
    }, 0)
  }))
  years <- generations / a$generations_per_year
  list(success = as.numeric(all(dose == 2L)),
       rp_recovery = 100 * (1 - mean(donor) / 2),
       donor_cM = sum(donor) / 2 * (sum(g$map$lengths) / length(donor)),
       linkage_drag = if (any(!is.na(drag))) mean(drag, na.rm = TRUE) else 0,
       generations = generations, years = years,
       plants = counter$plants, datapoints = counter$datapoints,
       cost = counter$plants * a$plant_cost + counter$datapoints * a$datapoint_cost)
}

simulate_backcross <- function(a, targets) {
  g <- bc_genome(a, targets)
  L <- g$map$n
  rp <- founder_pop(1L, L)
  donor <- founder_pop(2L, L)
  f1 <- mate(rp, 1L, donor, 1L, g$r)
  counter <- bc_counter()
  counter$plants <- 1
  all_t <- seq_along(g$targets)
  if (a$strategy == "simultaneous" || length(all_t) == 1) {
    best <- bc_stream(f1, rp, g, all_t, a, counter)
    final <- bc_fix(best, g, all_t, a$fix_plants, counter)
    return(bc_outcomes(final, g, a, counter, generations = 1 + a$bc_generations + 1))
  }
  # Parallel streams: convert and fix each target separately (at the same time).
  conversions <- lapply(all_t, function(i) {
    list(plant = bc_fix(bc_stream(f1, rp, g, i, a, counter), g, i, a$fix_plants, counter), need = i)
  })
  levels <- 0
  while (length(conversions) > 1) {
    levels <- levels + 1
    nxt <- list()
    for (j in seq(1, length(conversions), by = 2)) {
      if (j == length(conversions)) {
        nxt[[length(nxt) + 1]] <- conversions[[j]]
        next
      }
      p <- conversions[[j]]
      q <- conversions[[j + 1]]
      need <- c(p$need, q$need)
      n <- a$pyramid_plants
      pop <- list(h1 = gamete(p$plant, rep(1L, n), g$r), h2 = gamete(q$plant, rep(1L, n), g$r))
      counter$plants <- counter$plants + n
      keep <- bc_select(pop, g, need, 1L, counter, use_flanks = FALSE)
      nxt[[length(nxt) + 1]] <- list(plant = subset_pop(pop, keep), need = need)
    }
    conversions <- nxt
  }
  final <- bc_fix(conversions[[1]]$plant, g, all_t, a$fix_plants, counter)
  bc_outcomes(final, g, a, counter, generations = 1 + a$bc_generations + 1 + levels + 1)
}

# Theory -------------------------------------------------------------------------------------

#' Backcrossing by the numbers
#'
#' Classic expectations to check simulations against and to plan with.
#' `bc_expected_recovery()` gives the expected recurrent parent genome share
#' after `bc` backcross generations without background selection, 1 -
#' 0.5^(bc + 1). `min_population_size()` gives the smallest number of plants
#' that contains at least one plant of a genotype with frequency `p`, with
#' probability `prob`: log(1 - prob) / log(1 - p).
#'
#' @param bc Backcross generations (BC1 = 1).
#' @param p Frequency of the wanted genotype among progeny.
#' @param prob Wanted probability of finding at least one.
#' @return A numeric vector.
#' @examples
#' bc_expected_recovery(1:6)
#' # Plants needed to find a BC1 carrier recombinant on both sides of a gene
#' # with flanking markers 10 cM away (Haldane):
#' r <- 0.5 * (1 - exp(-2 * 10 / 100))
#' min_population_size(0.5 * r^2)
#' @export
bc_expected_recovery <- function(bc) 1 - 0.5^(bc + 1)

#' @rdname bc_expected_recovery
#' @export
min_population_size <- function(p, prob = 0.99) ceiling(log(1 - prob) / log(1 - p))
