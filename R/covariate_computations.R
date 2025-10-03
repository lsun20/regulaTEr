#' Preliminary covariate computations
#'
#' Performs preliminary computations only dependent on the covariates. Computes demeaned covariates and interaction terms baese on the specified treatment effect parameter.
#'
#' @param x A numeric vector; realizations of the treatment variable.
#' @param Z_uc A numeric matrix; realizations of the control covariates.
#' @param parameter A character string; Specifies which treatment effect ("ATE", "ATU", or "ATT") to compute. Defaults to "ATE".
#'
#' @returns A list containing the following components:
#' \itemize{
#'   \item \code{parameter}: A character string; the specified treatment effect ("ATE", "ATT", or "ATU").
#'   \item \code{n}: An integer; number of observations (rows).
#'   \item \code{x}: A numeric vector;realizations of the treatment variable.
#'   \item \code{Z_uc}: A numeric matrix; realizations of the control covariates.
#'   \item \code{Z_tilde}: A numeric matrix; the demeaned version of \code{Z_uc} according to the specified \code{parameter}.
#'   \item \code{Z_w}: A numeric matrix; realizations cof the interacted covariates.
#' }
covariate_comps <- function(x, Z_uc, parameter = "ATE") {
  n <- length(x)

  is_const_col <- function(col) all(col == 1)
  if(!any(apply(Z_uc, 2, is_const_col))) Z_uc <- cbind(Z_uc, 1)
  const_col <- apply(Z_uc, 2, is_const_col)

  const_col <- apply(Z_uc, 2, function(col) all(col == 1))
  w <- switch(parameter,
              ATT = x,
              ATU = 1 - x,
              rep(1, length(x)))  # default is ATE

  Z_tilde <- apply(Z_uc[, !const_col, drop = FALSE], 2, function(z) z - mean(w * z) / mean(w))
  Z_w <- x * Z_tilde

  list(parameter = parameter,
       n = n, x = x, Z_uc = Z_uc, Z_tilde = Z_tilde, Z_w = Z_w
  )
}
