#' Multidata Class
#'
#' The `multidata` class is an S7 class designed to represent multidimensional time series models, including the main time series and additional series (e.g., error standard deviations or related variables).
#'
#' @param times A numeric vector representing the time points of the time series.
#' @param series A numeric vector or matrix representing the main time series.
#' @param series_esd A numeric vector or matrix representing the additional series, such as error standard deviations or other related data.
#'
#' @section Validation Rules:
#' - `@times` and `@series` must be vectors or matrices (but not arrays).
#' - The lengths of `@times` and `@series` must be the same.
#' - If `@series_esd` is provided, it must be a vector or matrix with the same length as `@series`.
#'
#' @examples
#' # Create a multidata object
#' multidata_instance <- multidata(
#'   times = c(1, 2, 3, 4),
#'   series = c(10, 20, 15, 25),
#'   series_esd = c(1, 1.5, 1.2, 1.8)
#' )
#'
#' @export
multidata <- S7::new_class(
  "multidata",
  package = "iAR",
  properties = list(
    times = S7::class_numeric,
    series = S7::class_numeric,
    series_esd = S7::class_numeric
  )
)