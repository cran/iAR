#' `CiAR` Class
#'
#' Represents a complex irregular autoregressive (CiAR) time series model.
#' This class extends the `unidata` class and provides additional properties
#' for modeling, forecasting, and interpolation of complex-valued time series data.
#'
#' @param times A numeric vector representing the time points.
#' @param series A complex vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#' @param fitted_values A numeric vector containing the fitted values from the model.
#' @param kalmanlik A numeric value representing the Kalman likelihood of the model.
#' @param coef A numeric vector containing the estimated coefficients of the model.
#' @param tAhead A numeric value specifying the forecast horizon (default: 1).
#' @param forecast A numeric vector containing the forecasted values.
#' @param interpolated_values A numeric vector containing the interpolated values.
#' @param interpolated_times A numeric vector containing the times of the interpolated data points.
#' @param interpolated_series A numeric vector containing the interpolated series.
#' @param zero_mean A logical value indicating if the model assumes a zero-mean process (default: TRUE).
#' @param standardized A logical value indicating if the model assumes a standardized process (default: TRUE).
#'
#' @details
#' The `CiAR` class is designed to handle irregularly observed, complex-valued
#' time series data using an autoregressive approach. It extends the `unidata`
#' class to include functionalities specific to complex-valued data.
#'
#' Key features of the `CiAR` class include:
#' - Support for complex-valued time series data.
#' - Forecasting and interpolation functionalities for irregular time points.
#' - Assumptions of zero-mean and standardized processes, configurable by the user.
#'
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#'
#' @examples
#' o=iAR::utilities()
#' o<-gentime(o, n=200, distribution = "expmixture", lambda1 = 130, lambda2 = 6.5,p1 = 0.15, p2 = 0.85)
#' times=o@times
#' my_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#'
#' # Access properties
#' my_CiAR@coef
#'
#' @export
CiAR <- S7::new_class(
  "CiAR",
  parent = unidata,
  package = "iAR",
  properties = list(fitted_values = S7::class_numeric,
                    kalmanlik = S7::class_numeric,
                    coef = S7::class_numeric,
                    tAhead = S7::new_property(class_numeric, default = 1),
                    forecast = S7::class_numeric,
                    interpolated_values = S7::class_numeric,
                    interpolated_times = S7::class_numeric,
                    interpolated_series = S7::class_numeric,
                    zero_mean = S7::new_property(class_logical, default = TRUE),
                    standardized = S7::new_property(class_logical, default = TRUE))
)