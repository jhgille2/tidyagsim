map <- sim_genetic_map(2, 100, 2)
r <- tidyagsim:::recombination_fractions(map$chr, map$pos)

test_that("genetic map and recombination fractions", {
  expect_equal(map$n, 102)
  expect_equal(r[1], 0.5)
  expect_equal(r[52], 0.5)
  expect_equal(r[2], 0.5 * (1 - exp(-0.04)), tolerance = 1e-12)
})

founders <- tidyagsim:::founder_pop(1:2, map$n)
donor_share <- function(pop) mean(c(pop$h1, pop$h2) == 2)
het <- function(pop) mean(pop$h1 != pop$h2)

test_that("biparental populations have the expected genome composition", {
  set.seed(1)
  f2 <- tidyagsim:::make_family(founders, 1L, 2L, 400, 0, 2, r)
  expect_lt(abs(donor_share(f2) - (0.5)), 0.03)
  expect_lt(abs(het(f2) - (0.5)), 0.03)
  ril <- tidyagsim:::make_family(founders, 1L, 2L, 400, 0, 8, r)
  expect_lt(het(ril), 0.02)
  bc1 <- tidyagsim:::make_family(founders, 1L, 2L, 400, 1, 1, r)
  expect_lt(abs(donor_share(bc1) - (0.25)), 0.03)
  bc2f3 <- tidyagsim:::make_family(founders, 1L, 2L, 400, 2, 3, r)
  expect_lt(abs(donor_share(bc2f3) - (0.125)), 0.03)
  expect_lt(abs(het(bc2f3) - (0.125 * 0.25 * 2)), 0.03)
})

test_that("recombination between adjacent markers matches the map", {
  set.seed(2)
  bc1 <- tidyagsim:::make_family(founders, 1L, 2L, 4000, 1, 1, r)
  # Donor allele in the BC1 gamete from the F1 is h1; count switches in a 10 cM interval.
  h <- bc1$h1
  rec <- mean(h[1, ] != h[6, ])
  expect_equal(rec, 0.5 * (1 - exp(-0.2)), tolerance = 0.02)
})

test_that("MAGIC lines draw equally from each founder", {
  set.seed(3)
  f4 <- tidyagsim:::founder_pop(1:4, map$n)
  magic <- tidyagsim:::make_magic(f4, 4, 400, 6, r)
  share <- table(factor(c(magic$h1, magic$h2), levels = 1:4)) / length(magic$h1) / 2
  expect_equal(as.numeric(share), rep(0.25, 4), tolerance = 0.03)
  f8 <- tidyagsim:::founder_pop(1:8, map$n)
  magic8 <- tidyagsim:::make_magic(f8, 8, 200, 2, r)
  expect_equal(length(unique(c(magic8$h1))), 8)
})

test_that("marker scan LOD matches lm()", {
  set.seed(4)
  n <- 120
  markers <- matrix(sample(c(0L, 2L), 5 * n, replace = TRUE), 5)
  y <- 0.8 * markers[3, ] + rnorm(n)
  m5 <- list(chr = rep(1L, 5), pos = c(0, 10, 20, 30, 40))
  scan <- qtl_scan(markers, y, rep(1, n), m5, "marker", threshold = 3)
  lod_lm <- vapply(1:5, function(j) {
    (n / 2) * log10(sum(resid(lm(y ~ 1))^2) / sum(resid(lm(y ~ markers[j, ]))^2))
  }, numeric(1))
  expect_equal(scan$lod, lod_lm, tolerance = 1e-8)
  fam <- rep(1:3, each = 40)
  scan_f <- qtl_scan(markers, y, fam, m5, "family_marker", threshold = 3)
  lod_f <- vapply(1:5, function(j) {
    x <- markers[j, ]
    (n / 2) * log10(sum(resid(lm(y ~ factor(fam)))^2) / sum(resid(lm(y ~ factor(fam) + factor(fam):x))^2))
  }, numeric(1))
  expect_equal(scan_f$lod, lod_f, tolerance = 1e-8)
})

test_that("each population type simulates and reports outcomes", {
  for (type in c("F2", "RIL", "BCsFt", "NAM", "AB-NAM", "MAGIC")) {
    spec <- small_mapping(type = type, lines = 100, families = 4)
    out <- sim_simulate(spec, resolve_args(spec), 1L)[[1]]
    expect_true(out$power %in% c(0, 1), info = type)
    expect_equal(out$lines, 100, info = type)
    expect_true(out$threshold > 0, info = type)
    expect_true(out$cost > 0 && out$years >= 2, info = type)
  }
})

test_that("design checks reject impossible or over-budget designs", {
  check <- function(...) sim_check(small_mapping(...), resolve_args(small_mapping(...)))
  expect_null(check(type = "RIL", lines = 100))
  expect_match(check(type = "AB-NAM", backcrosses = 0), "backcross")
  expect_match(check(type = "NAM", lines = 20, families = 10), "too few")
  expect_match(check(type = "MAGIC", founders = "6"), "4 or 8")
  expect_match(check(type = "RIL", lines = 300, max_plots = 200), "plots exceed")
  expect_match(check(type = "RIL", generation = 1), "at least 2")
})

test_that("power rises with population size for a large QTL", {
  spec <- small_mapping(type = "RIL", lines = tune(), qtl_pve = 0.1, threshold = 3)
  res <- tune_grid(spec, sim_replicates(12), grid = data.frame(lines = c(40L, 300L)),
                   metrics = sim_metric_set(power))
  p <- collect_metrics(res)$mean
  expect_gt(p[2], p[1])
  expect_gt(p[2], 0.7)
})

test_that("designs share founders and QTL position on a seed", {
  a <- simulate_mapping_replicate(resolve_args(small_mapping(type = "RIL")), 7)
  b <- simulate_mapping_replicate(resolve_args(small_mapping(type = "NAM", families = 3)), 7)
  expect_equal(c(a$qtl_chr, a$qtl_pos), c(b$qtl_chr, b$qtl_pos))
})
