#' Estimate \code{C}.
#'
#' Estimate \code{C}, the heterogeneity bound, from the data. If the long regression includes
#' covariates with \code{NA} coefficients, use the trimmed data instead.
#'
#' @importFrom stats lm residuals
#' 
#' @param Y A numeric vector; realizations of the outcome variable.
#' @param x A numeric vector; realizations of the treatment variable.
#' @param Z_uc A numeric matrix; realizations of the control covariates.
#' @param Z_w A numeric matrix; realizations of the interacted covariates.
#' @param B A numeric matrix; covariance matrix of the demeaned \code{Z_uc}.
#' @param trimmed_data A list; output from \code{trim_data()}, used if the long regression has \code{NA} coefficients.
#' @param na_covariates A character vector; names of interaction covariates with \code{NA} coefficients with the
#' \code{Z_w} prefix removed obtained from running \code{find_long_NAs()}.
#'
#' @returns A list containing the square of the estimated \code{C} and its associated standard error, and CI.
#' 
cal_C <- function(Y, x, Z_uc, Z_w, B, trimmed_data, na_covariates) {
  if (length(na_covariates) > 0) {
    Y <- trimmed_data$Y
    x <- trimmed_data$x
    Z_uc <- trimmed_data$Z_uc
    Z_w <- trimmed_data$Z_w
    B <- trimmed_data$B
  }
  Z_combined <- cbind(Z_uc, x)

  proj <- lm(Z_w ~ Z_combined - 1)  # "-1" removes intercept to ensure exact projection
  Z_w_residualized <- residuals(proj)
  proj <- lm(Y ~ Z_combined - 1)  # "-1" removes intercept to ensure exact projection
  Y_residualized <- residuals(proj)
  model <- lm(Y_residualized ~ -1 + Z_w_residualized)  # "-1" removes intercept, dof will be different

  V <- diag(length(model$coefficients))
  bread <- qr.Q(model$qr) %*% backsolve(qr.R(model$qr),t(V), transpose = TRUE)

  B_matrix <- bread %*% B %*% t(bread )
  P_matrix <- Z_w_residualized %*% t(bread)

  M_matrix <- diag(length(Y)) - P_matrix
  diag_B <- diag(B_matrix)
  diag_M <- diag(M_matrix)
  diag_P <- diag(P_matrix)

  ratios <- diag_B / diag_M

  C_matrix_correction <- outer(ratios, ratios, '+')
  C_matrix <- B_matrix - 0.5*M_matrix * C_matrix_correction
  diag(C_matrix ) <- 0
  VAR_bc <- t(Y_residualized) %*% C_matrix %*% Y_residualized

  M_matrix_correction <- outer(diag_M, diag_M, '*') + diag_M^2
  C_matrix_double <- C_matrix / M_matrix_correction
  sigma2hats_resid <- Y_residualized * model$residuals
  VCATE_bc <- 2* t(sigma2hats_resid)%*%C_matrix_double^2%*%sigma2hats_resid
  C_sq <- as.numeric(VAR_bc)
  se <- sqrt(VCATE_bc)
  list(C_sq = C_sq, se = se, lb = C_sq -1.96*se, ub = C_sq + 1.96*se)
}
