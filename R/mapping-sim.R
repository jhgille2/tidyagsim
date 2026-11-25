# A fast simulator for QTL mapping populations ------------------------------------------------------
#
# Individuals are stored as two integer matrices (loci x individuals) holding the
# founder each allele descends from. All founders are inbred, so founder origin
# determines the marker genotype. Meiosis assumes no crossover interference
# (Haldane): the parity of crossovers in each interval between adjacent loci is an
# independent Bernoulli draw, which gives exact gametes at the simulated loci.

#' Simulate a genetic map
#'
#' Evenly spaced markers on chromosomes of equal genetic length. The defaults
#' (20 chromosomes of 115 cM, a marker per cM) approximate the size of
#' published soybean consensus linkage maps and a mid-density SNP panel.
#'
#' @param n_chromosomes Number of chromosomes.
#' @param chromosome_length Genetic length of each chromosome in cM (recycled).
#' @param marker_spacing Distance between markers in cM.
#' @return A list with `chr`, `pos` (cM) and `n`.
#' @examples
#' map <- sim_genetic_map(2, 50, 5)
#' str(map)
#' @export
sim_genetic_map <- function(n_chromosomes = 20, chromosome_length = 115, marker_spacing = 1) {
  lengths <- rep_len(chromosome_length, n_chromosomes)
  chr <- integer()
  pos <- numeric()
  for (k in seq_len(n_chromosomes)) {
    p <- seq(0, lengths[k], by = marker_spacing)
    chr <- c(chr, rep.int(k, length(p)))
    pos <- c(pos, p)
  }
  list(chr = chr, pos = pos, n = length(pos), lengths = lengths)
}

recombination_fractions <- function(chr, pos) {
  n <- length(pos)
  r <- rep(0.5, n)
  same <- c(FALSE, chr[-1] == chr[-n])
  d <- c(0, diff(pos))
  r[same] <- 0.5 * (1 - exp(-2 * d[same] / 100))
  r
}

founder_pop <- function(ids, L) {
  m <- matrix(rep(as.integer(ids), each = L), nrow = L)
  list(h1 = m, h2 = m)
}

gamete <- function(pop, which, r) {
  n <- length(which)
  h1 <- pop$h1[, which, drop = FALSE]
  h2 <- pop$h2[, which, drop = FALSE]
  # Every column starts at a chromosome boundary, where r = 0.5, so the parity of
  # the running count of switches is a fresh fair coin at each column start and no
  # per-column reset is needed.
  take2 <- (cumsum(stats::runif(length(r) * n) < r) %% 2L) == 1L
  h1[take2] <- h2[take2]
  h1
}

mate <- function(a, ia, b, ib, r) {
  list(h1 = gamete(a, ia, r), h2 = gamete(b, ib, r))
}

self_gen <- function(pop, r) {
  idx <- seq_len(ncol(pop$h1))
  mate(pop, idx, pop, idx, r)
}

# One biparental family: BCsFt with `s` backcrosses to founder p1 and filial generation t.
make_family <- function(founders, p1, p2, lines, s, t, r) {
  f1 <- mate(founders, p1, founders, p2, r)
  ones <- rep.int(1L, lines)
  if (s == 0) {
    pop <- mate(f1, ones, f1, ones, r)
    selfs <- t - 2
  } else {
    pop <- mate(f1, ones, founders, rep.int(p1, lines), r)
    for (k in seq_len(s - 1)) pop <- mate(pop, seq_len(lines), founders, rep.int(p1, lines), r)
    selfs <- t - 1
  }
  for (k in seq_len(max(selfs, 0))) pop <- self_gen(pop, r)
  pop
}

make_magic <- function(founders, n_founders, lines, t, r) {
  # Two-way F1s (one plant per pair, as seed is plentiful), then fresh plants at each
  # funnel level so every line traces to its own 4- or 8-way cross.
  pairs <- matrix(seq_len(n_founders), nrow = 2)
  level <- lapply(seq_len(ncol(pairs)), function(k) mate(founders, pairs[1, k], founders, pairs[2, k], r))
  ones <- rep.int(1L, lines)
  level <- lapply(level, function(p) list(h1 = p$h1[, ones, drop = FALSE], h2 = p$h2[, ones, drop = FALSE]))
  idx <- seq_len(lines)
  while (length(level) > 1) {
    level <- lapply(seq(1, length(level), by = 2), function(k) mate(level[[k]], idx, level[[k + 1]], idx, r))
  }
  pop <- level[[1]]
  for (k in seq_len(max(t - 1, 0))) pop <- self_gen(pop, r)
  pop
}

split_lines <- function(lines, families) {
  base <- lines %/% families
  base + (seq_len(families) <= lines %% families)
}

design_summary <- function(type, lines, families, s, t, founders) {
  if (type %in% c("MAGIC")) {
    funnel <- log2(founders)
    # 4-way: one 4-way cross per line; 8-way: two 4-way and one 8-way cross per line.
    crosses <- founders / 2 + lines * (founders / 2 - 1)
    generations <- funnel + t - 1
    advance <- lines * (funnel - 1 + t - 1)
  } else {
    fam <- if (type %in% c("NAM", "AB-NAM")) families else 1
    crosses <- fam * (1 + s)
    generations <- 1 + s + t - 1
    advance <- lines * (s + t - 1)
  }
  list(crosses = crosses, generations = generations, plant_generations = advance)
}

resolve_design <- function(a) {
  type <- a$type
  s <- switch(type, F2 = 0L, RIL = 0L, NAM = 0L, MAGIC = 0L, as.integer(a$backcrosses))
  t <- switch(type, F2 = 2L, as.integer(a$generation))
  families <- if (type %in% c("NAM", "AB-NAM")) as.integer(a$families) else 1L
  founders <- if (type == "MAGIC") as.integer(as.character(a$founders)) else if (families > 1) families + 1L else 2L
  list(type = type, s = s, t = t, families = families, founders = founders, lines = as.integer(a$lines))
}

#' Simulate one mapping population replicate
#'
#' The engine behind [mapping_population()]. Exposed for checking and teaching:
#' it returns the genotypes, phenotypes, QTL and scan as well as the outcomes.
#'
#' @param args A full argument list, as built by [mapping_population()].
#' @param seed Replicate seed.
#' @param keep_data Return the simulated genotypes, phenotypes and LOD profile.
#' @return A list of outcomes (`power` is 0 or 1 for one replicate), plus data
#'   when `keep_data = TRUE`.
#' @keywords internal
#' @export
simulate_mapping_replicate <- function(args, seed, keep_data = FALSE) {
  a <- args
  d <- resolve_design(a)
  set.seed(seed)
  map <- sim_genetic_map(a$n_chromosomes, a$chromosome_length, a$marker_spacing)
  # Design-independent draws first, so every design sees the same founders and QTL.
  max_founders <- 64L
  qchr <- sample.int(length(map$lengths), 1, prob = map$lengths)
  qpos <- stats::runif(1, 0, map$lengths[qchr])
  founder_alleles <- matrix(stats::rbinom(map$n * max_founders, 1, 0.5), nrow = map$n)
  qtl_bits <- stats::rbinom(max_founders, 1, a$qtl_donor_freq)
  bg_loci <- sort(sample.int(map$n, min(a$n_background, map$n)))
  bg_eff <- stats::rnorm(length(bg_loci))
  # Genome with the QTL inserted as an extra (untyped) locus.
  insert <- sum(map$chr < qchr | (map$chr == qchr & map$pos <= qpos))
  chr <- append(map$chr, qchr, after = insert)
  pos <- append(map$pos, qpos, after = insert)
  q_index <- insert + 1L
  marker_rows <- seq_along(chr)[-q_index]
  r <- recombination_fractions(chr, pos)
  L <- length(chr)
  if (d$founders > max_founders) abort(sprintf("At most %d founders are supported.", max_founders - 1))
  # QTL alleles: founder 1 (recurrent or common parent) carries the reference allele.
  qa <- qtl_bits[seq_len(d$founders)]
  qa[1] <- 0L
  if (d$type == "MAGIC") {
    qa[1] <- qtl_bits[1]
    if (length(unique(qa)) == 1) qa[1] <- 1L - qa[1]
  } else if (d$families == 1) {
    qa[2] <- 1L
  } else if (!any(qa[-1] == 1)) {
    qa[2] <- 1L
  }
  # Variance components on the plot scale (residual variance 1) in a reference biparental RIL.
  vg <- a$plot_h2 / (1 - a$plot_h2)
  vge <- a$gxe_ratio * vg
  bg_eff <- bg_eff * sqrt(2 * vg / sum(bg_eff^2))
  qtl_a <- sqrt(a$qtl_pve * (vg + vge + 1) / (1 - a$qtl_pve))
  # Build the population.
  founders <- founder_pop(seq_len(d$founders), L)
  if (d$type == "MAGIC") {
    pop <- make_magic(founders, d$founders, d$lines, d$t, r)
    family <- rep.int(1L, d$lines)
  } else {
    sizes <- split_lines(d$lines, d$families)
    parts <- lapply(seq_len(d$families), function(f) make_family(founders, 1L, f + 1L, sizes[f], d$s, d$t, r))
    pop <- list(h1 = do.call(cbind, lapply(parts, `[[`, "h1")), h2 = do.call(cbind, lapply(parts, `[[`, "h2")))
    family <- rep.int(seq_len(d$families), sizes)
  }
  n <- ncol(pop$h1)
  allele <- function(h) {
    rows <- row(h)
    out <- matrix(0L, nrow(h), ncol(h))
    is_marker <- rows != q_index
    marker_row <- match(rows, marker_rows)
    out[is_marker] <- founder_alleles[cbind(marker_row[is_marker], h[is_marker])]
    out[!is_marker] <- qa[h[!is_marker]]
    out
  }
  dosage <- allele(pop$h1) + allele(pop$h2)
  markers <- dosage[marker_rows, , drop = FALSE]
  g_bg <- as.vector(crossprod(markers[bg_loci, , drop = FALSE] - 1L, bg_eff))
  g_q <- qtl_a * (dosage[q_index, ] - 1)
  y <- g_bg + g_q + stats::rnorm(n, sd = sqrt(vge / a$locations + 1 / (a$locations * a$reps)))
  # Scan.
  model <- a$scan_model
  if (identical(model, "auto")) {
    model <- if (d$type == "MAGIC") "haplotype" else if (d$families > 1) "family_marker" else "marker"
  }
  scan <- qtl_scan(markers, y, family, map, model = model, n_perm = a$n_perm, alpha = a$alpha,
                   threshold = a$threshold, step = a$scan_step,
                   origin = if (model == "haplotype") list(h1 = pop$h1[marker_rows, , drop = FALSE],
                                                            h2 = pop$h2[marker_rows, , drop = FALSE],
                                                            founders = d$founders))
  call <- call_qtl(scan, qchr, qpos, a$window)
  info <- design_summary(d$type, d$lines, d$families, d$s, d$t, d$founders)
  plots <- n * a$locations * a$reps
  cost <- info$crosses * a$cross_cost + info$plant_generations * a$ssd_cost + n * a$genotype_cost +
    plots * a$plot_cost
  out <- list(power = call$detected, other_peaks = call$other_peaks, location_error = call$location_error,
              support_interval = call$support_interval, peak_lod = call$peak_lod, threshold = scan$threshold,
              realized_pve = stats::var(g_q) / stats::var(y), lines = n, plots = plots, genotyped = n,
              crosses = info$crosses, generations = info$generations,
              years = ceiling(info$generations / a$generations_per_year) + 1, cost = cost,
              qtl_chr = qchr, qtl_pos = qpos, segregating_markers = scan$segregating)
  if (keep_data) {
    out$data <- list(markers = markers, map = map, y = y, family = family, scan = scan, qtl_effect = qtl_a)
  }
  out
}

#' Genome scan for a single QTL
#'
#' Single-marker regression (`"marker"`), marker-within-family regression as in
#' nested association mapping joint linkage (`"family_marker"`), or regression
#' on founder-origin dosages (`"haplotype"`, for MAGIC). LOD scores compare
#' the residual sum of squares with and without the locus. The genome-wide
#' threshold is the `1 - alpha` quantile of the maximum LOD over `n_perm`
#' permutations of the phenotypes (within families when there are several).
#'
#' The haplotype model uses the true founder origin of each allele, so it is an
#' upper bound on what haplotype reconstruction (e.g. with qtl2) achieves.
#'
#' @param markers Marker dosage matrix (markers x lines, coded 0/1/2).
#' @param y Phenotypes.
#' @param family Family of each line.
#' @param map Genetic map from [sim_genetic_map()].
#' @param model Scan model.
#' @param n_perm Permutations for the threshold.
#' @param alpha Genome-wide type I error rate.
#' @param threshold A fixed LOD threshold instead of permutations.
#' @param step Scan every `step` cM instead of every marker (`NULL` for every marker).
#' @param origin Founder origin matrices for the haplotype model.
#' @return A list with the scan positions, `lod` and `threshold`.
#' @keywords internal
#' @export
qtl_scan <- function(markers, y, family, map, model = c("marker", "family_marker", "haplotype"), n_perm = 100,
                     alpha = 0.05, threshold = NULL, step = NULL, origin = NULL) {
  model <- match.arg(model)
  n <- length(y)
  keep <- seq_len(nrow(markers))
  if (!is.null(step) && step > 0) {
    grid_key <- paste(map$chr, floor(map$pos / step))
    keep <- keep[!duplicated(grid_key)]
  }
  perms <- if (is.null(threshold)) as.integer(n_perm) else 0L
  fam <- as.integer(factor(family))
  Y <- matrix(y, n, perms + 1)
  for (p in seq_len(perms)) {
    Y[, p + 1] <- if (max(fam) > 1) {
      unsplit(lapply(split(y, fam), function(v) v[sample.int(length(v))]), fam)
    } else y[sample.int(n)]
  }
  if (model == "marker" || (model == "family_marker" && max(fam) == 1)) {
    X <- t(markers[keep, , drop = FALSE])
    Xc <- sweep(X, 2, colMeans(X))
    sxx <- colSums(Xc^2)
    Yc <- sweep(Y, 2, colMeans(Y))
    sst <- colSums(Yc^2)
    sxy <- crossprod(Xc, Yc)
    ssr <- sxy^2 / sxx
    ssr[sxx < 1e-9, ] <- 0
  } else if (model == "family_marker") {
    ssr <- 0
    sst <- 0
    for (f in unique(fam)) {
      rows <- fam == f
      X <- t(markers[keep, rows, drop = FALSE])
      Xc <- sweep(X, 2, colMeans(X))
      sxx <- colSums(Xc^2)
      Yf <- Y[rows, , drop = FALSE]
      Yc <- sweep(Yf, 2, colMeans(Yf))
      sst <- sst + colSums(Yc^2)
      part <- crossprod(Xc, Yc)^2 / sxx
      part[sxx < 1e-9, ] <- 0
      ssr <- ssr + part
    }
  } else {
    if (is.null(origin)) abort("The haplotype model needs founder origins.")
    nf <- origin$founders
    Yc <- sweep(Y, 2, colMeans(Y))
    sst <- colSums(Yc^2)
    ssr <- matrix(0, length(keep), ncol(Y))
    for (j in seq_along(keep)) {
      m <- keep[j]
      D <- matrix(0, n, nf)
      D[cbind(seq_len(n), origin$h1[m, ])] <- 1
      D[cbind(seq_len(n), origin$h2[m, ])] <- D[cbind(seq_len(n), origin$h2[m, ])] + 1
      D <- sweep(D[, -nf, drop = FALSE], 2, colMeans(D[, -nf, drop = FALSE]))
      q <- qr(D)
      if (q$rank == 0) next
      Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
      ssr[j, ] <- colSums(crossprod(Q, Yc)^2)
    }
  }
  ssr <- pmin(ssr, rep(sst, each = nrow(ssr)) * (1 - 1e-12))
  lod <- (n / 2) * log10(rep(sst, each = nrow(ssr)) / (rep(sst, each = nrow(ssr)) - ssr))
  lod[!is.finite(lod)] <- 0
  thr <- if (is.null(threshold)) unname(stats::quantile(apply(lod[, -1, drop = FALSE], 2, max), 1 - alpha))
         else threshold
  list(chr = map$chr[keep], pos = map$pos[keep], lod = lod[, 1], threshold = thr,
       segregating = sum(apply(markers, 1, function(v) any(v != v[1]))))
}

call_qtl <- function(scan, qchr, qpos, window) {
  on_chr <- scan$chr == qchr
  near <- on_chr & abs(scan$pos - qpos) <= window
  peak_lod <- if (any(near)) max(scan$lod[near]) else 0
  detected <- as.numeric(peak_lod > scan$threshold)
  others <- !on_chr
  other_peaks <- as.numeric(any(others) && max(scan$lod[others]) > scan$threshold)
  location_error <- NA_real_
  support <- NA_real_
  if (detected == 1) {
    chr_lod <- scan$lod[on_chr]
    chr_pos <- scan$pos[on_chr]
    peak <- which.max(chr_lod)
    location_error <- abs(chr_pos[peak] - qpos)
    inside <- chr_lod >= chr_lod[peak] - 1.5
    left <- peak
    while (left > 1 && inside[left - 1]) left <- left - 1
    right <- peak
    while (right < length(chr_lod) && inside[right + 1]) right <- right + 1
    support <- chr_pos[right] - chr_pos[left]
  }
  list(detected = detected, other_peaks = other_peaks, location_error = location_error,
       support_interval = support, peak_lod = peak_lod)
}
