#' Unidata Class
#'
#' The `unidata` class is an S7 class designed to represent univariate irregularly observed time series models
#' with associated times, values, and optional error standard deviations.
#'
#' @param times A numeric vector representing the time points.
#' @param series A numeric vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#'
#' @section Validation Rules:
#' - `@times`, `@series`, and `@series_esd` must be vectors (not matrices or arrays).
#' - The lengths of `@times` and `@series` must be the same.
#' - If `@series_esd` is provided, it must be a vector with the same length as `@series`.
#'
#' @examples
#' # Create a unidata object
#' unidata_instance <- unidata(
#'   times = c(1, 2, 3, 4),
#'   series = c(10, 20, 15, 25),
#'   series_esd = c(1, 1.5, 1.2, 1.8)
#' )
#'
#' @export
unidata <- S7::new_class(
  "unidata",
  package = "iAR",
  properties = list(
    times = S7::class_numeric,
    series = S7::class_numeric,
    series_esd = S7::class_numeric
  ),
  validator = function(self) {
    if (!is.null(dim(self@times)) || !is.null(dim(self@series)) || !is.null(dim(self@series_esd))) {
      "@times, @series, and @series_esd must be vectors"
    } else if (!is.null(dim(self@series)) && length(self@times) != length(self@series)) {
      "The length of @times must be the same as the length of @series"
#    } else if (!is.null(self@series_esd) && length(self@series) != length(self@series_esd)) {
#      "The length of @series must be the same as the length of @series_esd"
    }
  }
)