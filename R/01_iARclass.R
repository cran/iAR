#' `iAR` Class
#'
#' Represents a univariate irregular autoregressive (iAR) time series model. 
#' This class extends the "unidata" class and includes additional properties 
#' for modeling, forecasting, and interpolation.
#'
#' @param times A numeric vector representing the time points.
#' @param series A numeric vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#' @param family A character string indicating the distribution family of the model (default: "norm").
#' @param fitted_values A numeric vector containing the fitted values from the model.
#' @param loglik A numeric value representing the log-likelihood of the model.
#' @param kalmanlik A numeric value representing the Kalman likelihood of the model.
#' @param coef A numeric vector containing the estimated coefficients of the model.
#' @param df A numeric value representing the degrees of freedom (`t` distribution).
#' @param sigma A numeric value representing the scale parameter (`t` distribution).
#' @param mean A numeric value representing the estimated mean of the model (`gamma` parameter).
#' @param variance A numeric value representing the estimated variance of the model (`gamma` parameter).
#' @param tAhead A numeric value specifying the forecast horizon (default: 1).
#' @param forecast A numeric vector containing the forecasted values.
#' @param interpolated_values A numeric vector containing the interpolated values.
#' @param interpolated_times A numeric vector containing the times of the interpolated data points.
#' @param interpolated_series A numeric vector containing the interpolated series.
#' @param zero_mean A logical value indicating if the model assumes a zero-mean process (default: TRUE).
#' @param standardized A logical value indicating if the model assumes a standardized process (default: TRUE).
#' @param hessian A logical value indicating whether the Hessian matrix is computed during estimation (default: FALSE).
#' @param summary A list containing the summary of the model fit, including diagnostics and statistical results.
#'
#' @details
#' The `iAR` class is designed to handle irregularly observed time series data using an 
#' autoregressive approach. It extends the "unidata" class to include additional 
#' modeling and diagnostic capabilities. Key functionalities include forecasting, 
#' interpolation, and model fitting.
#'
#' The class also supports advanced modeling features, such as:
#' - Different distribution families for the data (e.g., Gaussian, `t`-distribution).
#' - Optional computation of the Hessian matrix for parameter estimation.
#' - Standardized or zero-mean process assumptions.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @examples
#' # Create an `iAR` object
#' o=iAR::utilities()
#' o<-gentime(o, n=200, distribution = "expmixture", lambda1 = 130, lambda2 = 6.5,p1 = 0.15, p2 = 0.85)
#' times=o@times
#' my_iAR <- iAR(family = "norm", times = times, coef = 0.9,hessian=TRUE)
#'
#' my_iAR@family
#' my_iAR@coef
#'
#' @export
iAR <- S7::new_class(
  "iAR",
  parent = unidata,
  package = "iAR",
  properties = list(family = S7::new_property(S7::class_character,
                                              default = "norm"),
                    fitted_values = S7::class_numeric,
                    loglik = S7::class_numeric,
                    kalmanlik = S7::class_numeric,
                    coef = S7::class_numeric,
                    df = S7::class_numeric, # t
                    sigma = S7::class_numeric, # t
                    mean = S7::class_numeric, # gamma
                    variance = S7::class_numeric, #gamma
                    tAhead = S7::new_property(class_numeric, default = 1),
                    forecast = S7::class_numeric,
                    interpolated_values = S7::class_numeric,
                    interpolated_times = S7::class_numeric,
                    interpolated_series = S7::class_numeric,
                    zero_mean = S7::new_property(class_logical, default = TRUE), # no gamma
                    standardized = S7::new_property(class_logical, default = TRUE), # no gamma # no t
                    hessian = S7::new_property(class_logical, default = FALSE),
                    summary = S7::class_list)
)