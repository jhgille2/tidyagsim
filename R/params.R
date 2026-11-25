# dials-style parameter objects ---------------------------------------------------------

#' Parameter objects for simulation studies
#'
#' dials parameter functions for arguments of [mapping_population()] and other
#' specs. Use them to change a tuned range, exactly as `penalty()` is used in
#' tidymodels: `update(params, lines = n_lines(c(100, 600)))`.
#'
#' @param range A two-element vector with the inclusive lower and upper bound.
#' @param values Allowed values of a qualitative parameter.
#' @param trans A `trans` object from the scales package, or `NULL`.
#' @return A dials `param` object.
#' @name sim_params
#' @examples
#' n_lines()
#' n_lines(c(100, 400))
#' population_type(c("RIL", "NAM"))
NULL

#' @rdname sim_params
#' @export
n_lines <- function(range = c(50L, 1000L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(n_lines = "Lines in the population"), finalize = NULL)
}

#' @rdname sim_params
#' @export
n_families <- function(range = c(2L, 25L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(n_families = "Families (NAM donors)"), finalize = NULL)
}

#' @rdname sim_params
#' @export
backcross_generations <- function(range = c(0L, 3L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(backcross_generations = "Backcross generations (s)"), finalize = NULL)
}

#' @rdname sim_params
#' @export
filial_generation <- function(range = c(2L, 8L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(filial_generation = "Filial generation (t)"), finalize = NULL)
}

#' @rdname sim_params
#' @export
population_type <- function(values = c("F2", "RIL", "BCsFt", "NAM", "AB-NAM", "MAGIC")) {
  dials::new_qual_param(type = "character", values = values,
                        label = c(population_type = "Mapping population type"), finalize = NULL)
}

#' @rdname sim_params
#' @export
magic_founders <- function(values = c("4", "8")) {
  dials::new_qual_param(type = "character", values = as.character(values),
                        label = c(magic_founders = "MAGIC founders"), finalize = NULL)
}

#' @rdname sim_params
#' @export
n_locations <- function(range = c(1L, 4L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(n_locations = "Field locations"), finalize = NULL)
}

#' @rdname sim_params
#' @export
n_reps <- function(range = c(1L, 3L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(n_reps = "Replicates per location"), finalize = NULL)
}

#' @rdname sim_params
#' @export
qtl_pve <- function(range = c(0.02, 0.15), trans = NULL) {
  dials::new_quant_param(type = "double", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(qtl_pve = "Variance explained by the QTL"), finalize = NULL)
}

#' @rdname sim_params
#' @export
plot_heritability <- function(range = c(0.1, 0.8), trans = NULL) {
  dials::new_quant_param(type = "double", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(plot_heritability = "Plot-basis heritability"), finalize = NULL)
}

#' @rdname sim_params
#' @export
generations_per_year <- function(range = c(1L, 3L), trans = NULL) {
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(generations_per_year = "Generations per year"), finalize = NULL)
}

#' A quantitative parameter for any argument
#'
#' A shortcut for [dials::new_quant_param()] used by [sim_function()] and
#' [breeding_program()], where arguments have no built-in parameter function.
#'
#' @param range Inclusive lower and upper bounds.
#' @param type `"integer"` or `"double"`.
#' @param label Display label.
#' @return A dials `param` object.
#' @examples
#' sim_param(c(1, 10))
#' @export
sim_param <- function(range, type = NULL, label = "value") {
  type <- type %||% if (all(range == round(range))) "integer" else "double"
  if (type == "integer") range <- as.integer(range)
  dials::new_quant_param(type = type, range = range, inclusive = c(TRUE, TRUE), trans = NULL,
                         label = c(value = label), finalize = NULL)
}

#' @rdname sim_param
#' @param values Allowed values.
#' @export
sim_param_values <- function(values, label = "value") {
  if (is.logical(values)) {
    return(dials::new_qual_param(type = "logical", values = values, label = c(value = label), finalize = NULL))
  }
  dials::new_qual_param(type = "character", values = as.character(values), label = c(value = label),
                        finalize = NULL)
}
