# Plots ----------------------------------------------------------------------------------------------

#' Plot simulation results
#'
#' * `type = "marginals"`: mean +/- one standard error of each metric against
#'   each numeric parameter (qualitative parameters are shown as colour).
#' * `type = "pareto"`: the first two metrics against each other, with the
#'   Pareto front highlighted.
#' * `type = "race"` (racing results): designs remaining after each replicate.
#'
#' @param object A `sim_results` object or a run [sim_set()].
#' @param type Plot type.
#' @param metric Metric names to show; defaults to all (marginals) or the first two (pareto).
#' @param ... Not used.
#' @return A ggplot.
#' @exportS3Method ggplot2::autoplot
autoplot.sim_results <- function(object, type = c("marginals", "pareto", "race"), metric = NULL, ...) {
  rlang::check_installed("ggplot2")
  type <- match.arg(type)
  if (type == "pareto") return(plot_pareto(object, metric))
  if (type == "race") return(plot_race(object))
  long <- collect_metrics.sim_results(object, by_scenario = FALSE)
  if (!is.null(metric)) long <- long[long$.metric %in% metric, , drop = FALSE]
  ids <- object$params$id
  is_qual <- vapply(object$params$object, inherits, logical(1), "qual_param")
  num <- ids[!is_qual]
  colour <- if (any(is_qual)) ids[is_qual][1] else NULL
  long$.metric <- factor(long$.metric, levels = names(object$metrics))
  if (!length(num)) {
    x <- ids[1]
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[x]], y = .data$mean)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$mean - .data$std_err, ymax = .data$mean + .data$std_err)) +
      ggplot2::facet_wrap(ggplot2::vars(.data$.metric), scales = "free_y") +
      ggplot2::labs(y = "Mean (+/- 1 SE)")
    return(p + ggplot2::theme_bw())
  }
  pieces <- lapply(num, function(id) {
    d <- long
    d$.parameter <- id
    d$.value <- as.numeric(d[[id]])
    d
  })
  d <- do.call(rbind, pieces)
  aes_col <- if (is.null(colour)) ggplot2::aes() else ggplot2::aes(colour = .data[[colour]])
  ggplot2::ggplot(d, ggplot2::aes(x = .data$.value, y = .data$mean)) +
    aes_col +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean - .data$std_err, ymax = .data$mean + .data$std_err),
                           width = 0) +
    ggplot2::facet_grid(ggplot2::vars(.data$.metric), ggplot2::vars(.data$.parameter), scales = "free") +
    ggplot2::labs(x = NULL, y = "Mean (+/- 1 SE)") +
    ggplot2::theme_bw()
}

plot_pareto <- function(object, metric) {
  metrics <- if (inherits(object, "sim_set")) object$result[[1]]$metrics else object$metrics
  metric <- metric %||% utils::head(names(metrics), 2)
  if (length(metric) != 2) abort("A Pareto plot needs two metrics.")
  wide <- if (inherits(object, "sim_set")) collect_metrics.sim_set(object, by_scenario = FALSE, type = "wide")
          else collect_metrics.sim_results(object, by_scenario = FALSE, type = "wide")
  front <- if (inherits(object, "sim_set")) {
    pts <- cbind(metric_sign(object$result[[1]]$metrics[[metric[1]]]) * wide[[metric[1]]],
                 metric_sign(object$result[[1]]$metrics[[metric[2]]]) * wide[[metric[2]]])
    vapply(seq_len(nrow(pts)), function(i) all(is.finite(pts[i, ])) && !any(vapply(seq_len(nrow(pts)), function(j)
      j != i && all(is.finite(pts[j, ])) && all(pts[j, ] >= pts[i, ]) && any(pts[j, ] > pts[i, ]), logical(1))),
      logical(1))
  } else {
    wide$.config %in% pareto_front(object, metric)$.config
  }
  wide$front <- ifelse(front, "Pareto front", "Dominated")
  colour <- if (inherits(object, "sim_set")) ggplot2::aes(colour = .data$wflow_id) else ggplot2::aes()
  ggplot2::ggplot(wide, ggplot2::aes(x = .data[[metric[1]]], y = .data[[metric[2]]])) +
    ggplot2::geom_line(data = wide[front, , drop = FALSE][order(wide[[metric[1]]][front]), , drop = FALSE],
                       colour = "grey50", linetype = 2) +
    ggplot2::geom_point(colour) +
    ggplot2::geom_point(data = wide[front, , drop = FALSE], shape = 21, size = 3.5, colour = "black") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[metric[2]]] - .data[[paste0(metric[2], "_se")]],
                                        ymax = .data[[metric[2]]] + .data[[paste0(metric[2], "_se")]]),
                           width = 0, alpha = 0.4) +
    ggplot2::theme_bw()
}

plot_race <- function(object) {
  if (!nrow(object$trace) || !"kept" %in% names(object$trace)) abort("No racing trace in these results.")
  ggplot2::ggplot(object$trace, ggplot2::aes(x = .data$replicates, y = .data$kept)) +
    ggplot2::geom_step() + ggplot2::geom_point() +
    ggplot2::labs(x = "Replicates simulated", y = "Designs still in the race") +
    ggplot2::theme_bw()
}

#' @rdname autoplot.sim_results
#' @exportS3Method ggplot2::autoplot
autoplot.sim_set <- function(object, type = c("pareto", "ranks"), metric = NULL, ...) {
  rlang::check_installed("ggplot2")
  type <- match.arg(type)
  if (type == "pareto") return(plot_pareto(object, metric))
  ranks <- rank_results(object, metric[1])
  m <- metric[1] %||% names(object$result[[1]]$metrics)[1]
  ggplot2::ggplot(ranks, ggplot2::aes(x = .data$rank, y = .data[[m]], colour = .data$wflow_id)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[m]] - .data[[paste0(m, "_se")]],
                                        ymax = .data[[m]] + .data[[paste0(m, "_se")]]), width = 0) +
    ggplot2::labs(x = "Design rank", colour = NULL) +
    ggplot2::theme_bw()
}
