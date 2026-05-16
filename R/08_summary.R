#' Summary of iAR package models
#'
#' This method provides a summary for objects created by the iAR package,
#' including `iAR`, `CiAR`, `BiAR`, and `utilities` objects. For fitted model
#' objects, the summary reports model-specific information such as the model
#' family, coefficients, standard errors, t-values, and p-values.
#' The output is formatted to display the relevant statistics based on the model class.
#'
#' @name summary
#'
#' @param object An object of class `iAR`, `CiAR`, `BiAR`, or `utilities`.
#'   For fitted model objects, this contains model-specific parameters such as
#'   coefficients (`coef`), standard errors (`stderr`), t-values (`tvalue`),
#'   p-values (`pvalue`), and other class-specific components.
#' @param ... Additional arguments (unused).
#'
#'
#' @export
S7::method(summary, iAR) <- function(object,diagnostic=TRUE,...) {
  x   <- object
  if(length(x@loglik) == 0 && length(x@kalmanlik) == 0) {
    stop("The model has not been fitted. Use loglik() or kalman() before summary().",call. = FALSE)
  }
  n   <- length(x@series)
  k   <- 1
  has_hess <- length(x@summary) > 0 && !is.null(x@summary$stderr)
  if (length(x@loglik)    > 0) nll <- x@loglik
  if (length(x@kalmanlik) > 0) nll <- 0.5 * n * (log(2*pi) + x@kalmanlik)
  AIC <- 2*nll + 2*k
  BIC <- 2*nll + log(n)*k
  cat("iAR model\n\n")
  if (x@family == "norm") {
    cat("Family:\n ", x@family, "\n\n")
    cat("Coefficients:\n")
    if (has_hess) {
      cat(paste("    ", "Estimate St. Error t value Pr(>|t|)\n"))
      cat("phi", "     ", format(round(x@coef,2),nsmall=2), "      ",
          .fmt_se(x@summary$stderr), "   ",
          .fmt_se(x@summary$tvalue), "     ",
          .fmt_se(x@summary$pvalue), "\n", sep = "")
    } else {
      cat("    Estimate\nphi     ", format(round(x@coef,2),nsmall=2), "\n", sep = "")
      cat("(Set hessian=TRUE before loglik()/kalman() to obtain standard errors)\n")
    }
  }
  if (x@family == "t") {
    k <- 2; AIC <- 2*nll + 2*k; BIC <- 2*nll + log(n)*k
    cat("Family:\n ", x@family, "\n\nd.f.:\n ", x@df, "\n\nCoefficients:\n")
    if (has_hess) {
      cat(paste("      ", "Estimate St. Error t value Pr(>|t|)\n"))
      cat("phi", "       ", format(round(x@coef,2),nsmall=2), "      ",
          .fmt_se(x@summary$stderr[1]), "   ", .fmt_se(x@summary$tvalue[1]), "     ",
          .fmt_se(x@summary$pvalue[1]), "\n", sep = "")
      cat("sigma", "     ", format(round(x@sigma,2),nsmall=2), "      ",
          .fmt_se(x@summary$stderr[2]), "    ", .fmt_se(x@summary$tvalue[2]), "     ",
          .fmt_se(x@summary$pvalue[2]), "\n", sep = "")
    } else {
      cat("    Estimate\nphi     ", format(round(x@coef,2),nsmall=2),
          "\nsigma   ", format(round(x@sigma,2),nsmall=2), "\n", sep = "")
      cat("(Set hessian=TRUE before loglik() to obtain standard errors)\n")
    }
  }
  if (x@family == "gamma") {
    k <- 3; AIC <- 2*nll + 2*k; BIC <- 2*nll + log(n)*k
    cat("Family:\n ", x@family, "\n\nCoefficients:\n")
    if (has_hess) {
      cat(paste("         ", "Estimate St. Error t value Pr(>|t|)\n"))
      cat("phi", "          ", format(round(x@coef,2),nsmall=2), "      ",
          .fmt_se(x@summary$stderr[1]), "   ", .fmt_se(x@summary$tvalue[1]), "     ",
          .fmt_se(x@summary$pvalue[1]), "\n", sep = "")
      cat("mean", "         ", format(round(x@mean,2),nsmall=2), "      ",
          .fmt_se(x@summary$stderr[2]), "   ", .fmt_se(x@summary$tvalue[2]), "     ",
          .fmt_se(x@summary$pvalue[2]), "\n", sep = "")
      cat("variance", "     ", format(round(x@variance,2),nsmall=2), "      ",
          .fmt_se(x@summary$stderr[3]), "    ", .fmt_se(x@summary$tvalue[3]), "     ",
          .fmt_se(x@summary$pvalue[3]), "\n", sep = "")
    } else {
      cat("    Estimate\nphi      ", format(round(x@coef,2),nsmall=2),
          "\nmean     ", format(round(x@mean,2),nsmall=2),
          "\nvariance ", format(round(x@variance,2),nsmall=2), "\n", sep = "")
      cat("(Set hessian=TRUE before loglik() to obtain standard errors)\n")
    }
  }
  cat("\nInformation criteria:\n")
  cat("AIC:", format(round(AIC,2), nsmall=2), "\n")
  cat("BIC:", format(round(BIC,2), nsmall=2), "\n")
  if(diagnostic)
  {
    x <- fit(x)
    residuals <- x@series - x@fitted_values
    cat("\n", sep = "")
    cat("Residual diagnostics:", "\n", sep = "")
    cat("ACF lag 1:",
      format(round(acf(residuals, plot = FALSE)$acf[2], 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=1):",
      format(round(Box.test(residuals, lag = 1, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=10):",
      format(round(Box.test(residuals, lag = 10, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
  }
  if (has_hess) {
    if (!isTRUE(tryCatch({ chol(x@summary$hess); TRUE }, error = function(e) FALSE)))
      warning("Hessian matrix is not positive definite or not invertible.", call. = FALSE)
    if (any(!is.finite(x@summary$stderr)))
      warning("Some standard errors are non-finite.", call. = FALSE)
    if (!is.null(x@summary$convergence) && x@summary$convergence != 0)
      warning(paste0("Optimization may not have converged (code = ", x@summary$convergence, ")."), call. = FALSE)
  }
}

S7::method(summary, CiAR) <- function(object,diagnostic=TRUE,...) {
  x <- object
  if(length(x@kalmanlik) == 0) {
    stop("The model has not been fitted. Use kalman() before summary().",call. = FALSE)}
  n <- length(x@series)
  k <- 2
  has_hess <- length(x@summary) > 0 && !is.null(x@summary$stderr)
  nll <- 0.5 * n * (log(2*pi) + x@kalmanlik)
  AIC <- 2*nll + 2*k
  BIC <- 2*nll + log(n)*k
  cat("CiAR model\n\nCoefficients:\n")
  if (has_hess) {
    cat(paste("    ", "Estimate St. Error t value Pr(>|t|)\n"))
    cat("phiR", "     ", format(round(x@coef[1],2),nsmall=2), "      ",
        .fmt_se(x@summary$stderr[1]), "   ", .fmt_se(x@summary$tvalue[1]), "     ",
        .fmt_se(x@summary$pvalue[1]), "\n", sep = "")
    cat("phiI", "     ", format(round(x@coef[2],2),nsmall=2), "      ",
        .fmt_se(x@summary$stderr[2]), "    ", .fmt_se(x@summary$tvalue[2]), "     ",
        .fmt_se(x@summary$pvalue[2]), "\n", sep = "")
  } else {
    cat("     Estimate\nphiR     ", format(round(x@coef[1],2),nsmall=2),
        "\nphiI     ", format(round(x@coef[2],2),nsmall=2), "\n", sep = "")
    cat("(Set hessian=TRUE before kalman() to obtain standard errors)\n")
  }
  cat("\nInformation criteria:\n")
  cat("AIC:", format(round(AIC,2), nsmall=2), "\n")
  cat("BIC:", format(round(BIC,2), nsmall=2), "\n")
  cat("\n", sep = "")
  if(diagnostic)
  {
    x <- fit(x)
    residuals <- x@series - x@fitted_values
    cat("Residual diagnostics:", "\n", sep = "")
    cat("ACF lag 1:",
      format(round(acf(residuals, plot = FALSE)$acf[2], 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=1):",
      format(round(Box.test(residuals, lag = 1, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=10):",
      format(round(Box.test(residuals, lag = 10, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
  }
  if (has_hess) {
    if (!isTRUE(tryCatch({ chol(x@summary$hess); TRUE }, error = function(e) FALSE)))
      warning("Hessian matrix is not positive definite or not invertible.", call. = FALSE)
    if (any(!is.finite(x@summary$stderr)))
      warning("Some standard errors are non-finite.", call. = FALSE)
    if (!is.null(x@summary$convergence) && x@summary$convergence != 0)
      warning(paste0("Optimization may not have converged (code = ", x@summary$convergence, ")."), call. = FALSE)
  }
}

S7::method(summary, BiAR) <- function(object,diagnostic = TRUE,...) {
  x <- object
  if(length(x@kalmanlik) == 0) {
    stop("The model has not been fitted. Use kalman() before summary().",call. = FALSE)}
  n <- dim(x@series)[1]
  k <- 2
  has_hess <- length(x@summary) > 0 && !is.null(x@summary$stderr)
  nll <- 0.5 * n * (log(2*pi) + x@kalmanlik)
  AIC <- 2*nll + 2*k
  BIC <- 2*nll + log(n)*k
  cat("BiAR model\n\nCoefficients:\n")
  if (has_hess) {
    cat(paste("    ", "Estimate St. Error t value Pr(>|t|)\n"))
    cat("phiR", "     ", format(round(x@coef[1],2),nsmall=2), "      ",
        .fmt_se(x@summary$stderr[1]), "   ", .fmt_se(x@summary$tvalue[1]), "     ",
        .fmt_se(x@summary$pvalue[1]), "\n", sep = "")
    cat("phiI", "     ", format(round(x@coef[2],2),nsmall=2), "      ",
        .fmt_se(x@summary$stderr[2]), "    ", .fmt_se(x@summary$tvalue[2]), "     ",
        .fmt_se(x@summary$pvalue[2]), "\n", sep = "")
  } else {
    cat("     Estimate\nphiR     ", format(round(x@coef[1],2),nsmall=2),
        "\nphiI     ", format(round(x@coef[2],2),nsmall=2), "\n", sep = "")
    cat("(Set hessian=TRUE before kalman() to obtain standard errors)\n")
  }
  cat("\nInformation criteria:\n")
  cat("AIC:", format(round(AIC,2), nsmall=2), "\n")
  cat("BIC:", format(round(BIC,2), nsmall=2), "\n\n")
  cat("\n", sep = "")
  if(diagnostic)
  {
    x <- fit(x)
    residuals <- x@series - x@fitted_values
    cat("Residual diagnostics:", "\n", sep = "")
    cat("ACF lag 1 residual 1:",
      format(round(acf(residuals[,1], plot = FALSE)$acf[2], 2), nsmall = 2),
      "\n", sep = " ")
    cat("ACF lag 1 residual 2:",
      format(round(acf(residuals[,2], plot = FALSE)$acf[2], 2), nsmall = 2),
      "\n", sep = " ")
    ccf_res <- ccf(residuals[,1], residuals[,2], plot = FALSE, lag.max = 1)
    cat("CCF between residual 1(t) and residual 2(t-1):",
      format(round(ccf_res$acf[1], 2), nsmall = 2), "\n", sep = " ")
    cat("CCF between residual 1(t-1) and residual 2(t):",
      format(round(ccf_res$acf[3], 2), nsmall = 2), "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=1) residual 1:",
      format(round(Box.test(residuals[,1], lag = 1, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=1) residual 2:",
      format(round(Box.test(residuals[,2], lag = 1, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=10) residual 1:",
      format(round(Box.test(residuals[,1], lag = 10, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
    cat("Ljung-Box test p-value (lag=10) residual 2:",
      format(round(Box.test(residuals[,2], lag = 10, type = "Ljung-Box")$p.value, 2), nsmall = 2),
      "\n", sep = " ")
  }
  if (has_hess) {
    if (!isTRUE(tryCatch({ chol(x@summary$hess); TRUE }, error = function(e) FALSE)))
      warning("Hessian matrix is not positive definite or not invertible.", call. = FALSE)
    if (any(!is.finite(x@summary$stderr)))
      warning("Some standard errors are non-finite.", call. = FALSE)
    if (!is.null(x@summary$convergence) && x@summary$convergence != 0)
      warning(paste0("Optimization may not have converged (code = ", x@summary$convergence, ")."), call. = FALSE)
  }
}

S7::method(summary, utilities) <- function(object,...) {
  print(object@summary)}
