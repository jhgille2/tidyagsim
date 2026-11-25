# Replicates and scenarios: the rsample analogue -----------------------------------------

#' Paired replicate seeds
#'
#' In a simulation study the unit of replication is a random seed. Every
#' candidate design is run on the *same* seeds (common random numbers), so
#' differences between designs are paired: the same founders, QTL positions and
#' field errors are drawn for each design wherever the designs allow it.
#' `sim_replicates()` plays the role of `rsample::vfold_cv()`.
#'
#' @param times Number of replicates.
#' @param seed First seed; replicate `i` uses `seed + i - 1`.
#' @param scenarios Optional assumption scenarios from [sim_scenarios()]. Each
#'   design is run on every replicate in every scenario.
#' @return A tibble of class `sim_rset` with columns `id` and `seed`, and the
#'   scenarios as an attribute.
#' @examples
#' sim_replicates(20)
#' sim_replicates(10, scenarios = sim_scenarios(qtl_pve = c(0.03, 0.08)))
#' @export
sim_replicates <- function(times = 10, seed = 1, scenarios = NULL) {
  if (!is.numeric(times) || length(times) != 1 || times < 1 || times != round(times)) {
    abort("`times` must be a positive whole number.")
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed < 1 || seed + times - 1 > .Machine$integer.max) {
    abort("`seed` must be a positive whole number.")
  }
  width <- nchar(as.character(times))
  out <- tibble::tibble(id = sprintf(paste0("Rep%0", width, "d"), seq_len(times)),
                        seed = as.integer(seed + seq_len(times) - 1))
  scenarios <- scenarios %||% sim_scenarios()
  if (!inherits(scenarios, "sim_scenarios")) abort("`scenarios` must come from sim_scenarios().")
  structure(out, scenarios = scenarios, class = c("sim_rset", class(out)))
}

#' Assumption scenarios
#'
#' Crosses assumption factors into named scenarios, for things a breeder cannot
#' choose but must plan around (heritability, marker cost, QTL effect size).
#' Each name must be an argument of the spec being tuned; the scenario value
#' overrides the spec's fixed value.
#'
#' @param ... Named vectors. Scalars are held fixed but still recorded.
#' @return A tibble of class `sim_scenarios` with a `.scenario` label column and
#'   one column per factor.
#' @examples
#' sim_scenarios(plot_h2 = c(0.3, 0.6), genotype_cost = c(15, 30))
#' @export
sim_scenarios <- function(...) {
  factors <- list(...)
  if (!length(factors)) {
    out <- tibble::tibble(.scenario = "base")
    return(structure(out, class = c("sim_scenarios", class(out))))
  }
  if (is.null(names(factors)) || any(names(factors) == "")) abort("Name every scenario factor.")
  grid <- expand.grid(factors, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  labels <- vapply(seq_len(nrow(grid)), function(i) {
    paste(names(grid), vapply(grid[i, , drop = FALSE], function(v) as.character(v[[1]]), ""), sep = "=",
          collapse = ", ")
  }, character(1))
  out <- tibble::as_tibble(cbind(data.frame(.scenario = labels, stringsAsFactors = FALSE), grid))
  if (anyDuplicated(out$.scenario)) abort("Scenario factors produce duplicate labels.")
  structure(out, class = c("sim_scenarios", class(out)))
}

scenario_list <- function(rset) {
  sc <- attr(rset, "scenarios") %||% sim_scenarios()
  lapply(seq_len(nrow(sc)), function(i) {
    values <- as.list(sc[i, setdiff(names(sc), ".scenario"), drop = FALSE])
    list(id = sc$.scenario[i], values = values)
  })
}

as_rset <- function(resamples) {
  if (is.null(resamples)) return(sim_replicates(10))
  if (is.numeric(resamples) && length(resamples) == 1) return(sim_replicates(resamples))
  if (!inherits(resamples, "sim_rset")) abort("`resamples` must come from sim_replicates() (or be a number).")
  resamples
}

#' @export
print.sim_rset <- function(x, ...) {
  sc <- attr(x, "scenarios")
  cat(sprintf("# Simulation replicates: %d paired seeds (%d..%d) x %d scenario(s)\n",
              nrow(x), min(x$seed), max(x$seed), nrow(sc)))
  if (nrow(sc) > 1) cat("# Scenarios:", paste(sc$.scenario, collapse = " | "), "\n")
  print(tibble::as_tibble(unclass_rset(x)), ...)
  invisible(x)
}

unclass_rset <- function(x) {
  attr(x, "scenarios") <- NULL
  class(x) <- setdiff(class(x), "sim_rset")
  x
}

#' @export
`[.sim_rset` <- function(x, i, j, ...) {
  sc <- attr(x, "scenarios")
  out <- NextMethod()
  if (missing(j) && all(c("id", "seed") %in% names(out))) {
    attr(out, "scenarios") <- sc
    class(out) <- unique(c("sim_rset", class(out)))
  }
  out
}
