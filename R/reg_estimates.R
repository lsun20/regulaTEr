#' Matrix square root
#'
#' Computes the (approximate) principal square root \eqn{B} of a matrix \eqn{A},
#' so that \eqn{BB \approx A}, via Schur decomposition.
#'
#' @param A A numeric matrix.
#'
#' @return A numeric matrix; representing the (approximate) principal square root of \code{A}.
sqrt_mat <- function(A) {
  schur_decomp <- eigen(A)
  U <- schur_decomp$vectors
  T <- schur_decomp$values
  sqrt_T <- diag(sqrt(T))
  return(U %*% sqrt_T %*% solve(U))
}

#' \code{NA} coefficients in long regression detector
#'
#' Returns the names of the \code{NA} coefficients in the interaction covariates of a given long regression: these
#' are the coefficients contributing to lack of overlap and hence a long regression with \code{NA} coefficients.
#' 
#' @importFrom stats coef
#' 
#' @param long_reg An \code{lm} object; the long regression object, \code{lm(Y ~ -1 + x + Z_conf + Z_w)}.
#'
#' @returns A character vector; containing the names of interaction covariates with \code{NA} coefficients with the
#' \code{Z_w} prefix removed. If there are no such covariates, returns an empty character vector.
#' 
find_long_NAs <- function(long_reg) {
  coefficients <- coef(long_reg)
  covar_to_remove <- character()

  for (name in names(coefficients)) {
    if (is.na(coefficients[name]) && startsWith(name, "Z_w")) {
      clean_name <- sub("^Z_w", "", name)
      covar_to_remove <- c(covar_to_remove, clean_name)
    }
  }
  return (covar_to_remove)
}

#' Ridge long regression
#'
#' Runs a ridge long regression with a penalty determined by 10-fold cross-validation and computes
#' the residuals and effective degrees of freedom. Used when the long regression has \code{NA} coefficients.
#'
#' @param Y A numeric vector; realizations of the outcome variable.
#' @param x A numeric vector; realizations of the treatment variable.
#' @param Z_conf A numeric matrix; realizations of the confounder covariates.
#' @param Z_w A numeric matrix; realizations cof the interacted covariates.
#' @param B A numeric matrix; covariance matrix of the demeaned \code{Z_uc}.
#'
#' @returns A list containing the following components:
#' \itemize{
#'   \item \code{ridge}: An \code{glmnet} object; fitted ridge regression model, an object of class \code{glmnet}.
#'   \item \code{residuals}: A numeric vector; residuals.
#'   \item \code{df}: A numeric vector; the effective degrees of freedom for the ridge regression model,
#'               computed as \code{n - sum(d^2 / (d^2 + lambda))}, where \code{d} are the singular
#'               values of the design matrix \code{X} and \code{lambda} is the optimal penalty
#'               determined by cross-validation.
#' }
run_long_ridge <- function(Y, x, Z_conf, Z_w, B) {
  Z_w_S <- Z_w %*% solve(sqrt_mat(B))
  X <- cbind(x, Z_conf, Z_w_S)

  cv <- glmnet::cv.glmnet(X, Y, alpha = 0, nfolds = 10)

  ridge <- glmnet::glmnet(X, Y, alpha = 0, lambda = cv$lambda.min)
  residuals <- (Y - predict(ridge, newx = X))

  d <- svd(X)$d
  df <- length(Y) - (sum(d^2 / (d^2 + cv$lambda.min)))
  list(ridge = ridge, residuals = residuals, df = df)
}

#' Long or ridge regression residuals
#'
#' Returns residuals, either from the long regression or ridge long regression,
#' depending on whether the long regression contains any \code{NA} coefficients.
#'
#' @param long_NAs A character vector; names of interaction covariates with \code{NA} coefficients, with the \code{Z_w} prefix removed. Possibly empty.
#' @param long_reg An \code{lm} object; the long regression object, \code{lm(Y ~ -1 + x + Z_uc + Z_w)}.
#' @param ridgereg A \code{glmnet} object; the ridge regression object obtained from running \code{run_long_ridge()}.
#'
#' @returns A list containing the following components:
#' \itemize{
#'   \item \code{residuals}: A numeric vector; residuals, either from the long regression or ridge long regression.
#'   \item \code{df}: A numeric vector; degrees of freedom corresponding to the long regression or ridge long regression.
#' }
determine_residuals <- function(long_NAs, long_reg, ridgereg) {
  if (length(long_NAs) > 0) {
    res <- ridgereg$residuals
    df <- ridgereg$df
  } else {
    res <- long_reg$residuals
    df <- long_reg$df
  }
  list (res = res, df = df)
}

#' Short, long, and ridge regressions
#'
#' Estimates long and short regressions, as well as a trimmed long and ridge regression
#' (via \code{run_long_ridge()}) if \code{NA} coefficients exist.
#'
#' @param Y A numeric vector; realizations of the outcome variable.
#' @param x A numeric vector; realizations of the treatment variable.
#' @param covariate_comps A list; covariate-related components returned by \code{covariate_comps()}.
#' @param Z_conf A numeric matrix; realizations of the confounder covariates. Defaults to \code{Z_uc} if \code{NULL}.
#' @param trimmed_data A list; contains the matrices/vectors for the outcome (\code{Y}), 
#' treatment indicator (\code{x}), and interaction variables (\code{Z_w}), plus either 
#' covariates (\code{Z_uc}) or confounders (\code{Z_conf}). Used when trimming is necessary due to lack of 
#' overlap and the user wants to specify the trimmed data directly. Defaults to \code{NULL}, 
#' in which case the function trims the data automatically.
#' 
#' @returns A list with the following components:
#' \itemize{
#'   \item \code{multi_results}: A list; contains (if applicable):
#'     \itemize{
#'       \item \code{ridgereg}: A \code{glmnet} object; the ridge regression object (or \code{NULL} if unused).
#'       \item \code{trim_reg}: An \code{lm} object; the trimmed long regression object.
#'       \item \code{res_long_prop_trim}: A numeric vector; residuals from projecting \code{x} on trimmed covariates.
#'       \item \code{long_wghts_trim}: A numeric vector; weights from trimmed long regression.
#'       \item \code{long_wghts_trim_full}: A numeric vector; weights from trimmed long regression, but with zeros at trimmed rows.
#'       \item \code{rows_to_delete}: An integer vector; trimmed row indices.
#'     }
#'   \item \code{covariate_comps}: A list; the original, unchanged input to the function.
#'   \item \code{long_reg}: A \code{lm} object; \code{Y ~ - 1 + x + Z_conf + Z_w}.
#'   \item \code{long_NAs}: An integer vector; indices for coefficients dropped in the long regression.
#'   \item \code{num_covar_long}: An integer; number of covariates in the long regression.
#'   \item \code{num_covar_trim}: An integer; number of covariates in the trimmed regression.
#'   \item \code{shortreg}: An \code{lm} object; \code{Y ~ - 1 + x + Z_conf}.
#'   \item \code{residuals_df}: A list; obtained from \code{determine_residuals()} with:
#'     \itemize{
#'       \item \code{res}: A numeric vector; residual, either from the long or ridge long regression.
#'       \item \code{df}: An integer; residual degrees of freedom
#'     }
#'   \item \code{sigma}: A numeric vector; estimated residual standard deviation.
#'   \item \code{wghts}: A numeric matrix; the matrix cross-product of \code{short_wghts} and \code{Z_w}.
#'   \item \code{res_long_prop}: A numeric vector; residuals from projecting \code{x} onto \code{Z_conf} and \code{Z_w}.
#'   \item \code{long_wghts}: A numeric vector; weights used in the long regression specification.
#'   \item \code{num_covar_long}: An integer; number of covariates in the long regression.
#'   \item \code{res_prop}: A numeric vector; residuals from projecting \code{x} onto \code{Z_conf}.
#'   \item \code{short_wghts}: A numeric vector; weights used in the short regression specification.
#'   \item \code{short_prop}: An \code{lm} object; corresponds to the regression of \code{x} on \code{Z_conf}.
#'   \item \code{resid_Z_w_Z_uc}: A numeric matrix; residualized \code{Z_w} after projection onto \code{Z_conf}.
#'   \item \code{Z_w_resid_Z_w_Z_uc}: A numeric matrix; matrix cross-product of residualized \code{Z_w}.
#'   \item \code{Z_w_res_prop}: A numeric vector; matrix cross-product of residualized \code{Z_w} with \code{res_prop}.
#'   \item \code{B}: A numeric matrix; the covariance matrix of the demeaned \code{Z_uc}.
#'   \item \code{inv_B}: A numeric matrix; the inverse of \code{B}.
#' }
reg_estimates <- function(Y, x, covariate_comps, Z_conf = NULL, trimmed_data = NULL) {
  with(covariate_comps, {
    if (is.null(Z_conf)) {
      Z_conf <- Z_uc
    }
    
    B <- crossprod(Z_tilde) / n

    res_long_prop <- residuals(lm(x ~ -1 + Z_conf + Z_w))
    long_wghts <- res_long_prop / (crossprod(x, res_long_prop)[1, 1])
    
    short_prop <- lm(x ~ -1 + Z_conf)
    res_prop <- residuals(short_prop)
    short_wghts <- res_prop / (crossprod(x, res_prop)[1, 1])

    resid_Z_w_Z_uc <- residuals(lm(Z_w ~ -1 + Z_conf))
    Z_w_resid_Z_w_Z_uc <- crossprod(Z_w, resid_Z_w_Z_uc)
    Z_w_res_prop <- crossprod(Z_w, res_prop)

    long_reg <- lm(Y ~ -1 + x + Z_conf + Z_w)
    num_covar_long <- ncol(Z_conf) + ncol(Z_w)
    long_NAs <- find_long_NAs(long_reg)

    multi_comps <- function() {
      ridgereg <- run_long_ridge(Y, x, Z_conf, Z_w, B)
      if (is.null(trimmed_data)) {
        if(!isTRUE(all.equal(Z_conf, Z_uc))) stop("The function `trim_data()` can only be used when `Z_conf = Z_uc`. User should provide `trimmed_data` instead.")
        trimmed_data <- trim_data(cbind(Y, x, Z_uc), long_NAs, parameter)
      }
      

      if (is.null(trimmed_data$Z_conf)) trimmed_data$Z_conf <- trimmed_data$Z_uc
      trimmed_data$Z_w <- covariate_comps(trimmed_data$x, trimmed_data$Z_uc, parameter = parameter)$Z_w

      num_covar_trim <- ncol(trimmed_data$Z_conf) + ncol(trimmed_data$Z_w)
      trim_reg <- lm(trimmed_data$Y ~ -1 + trimmed_data$x + trimmed_data$Z_conf + trimmed_data$Z_w)
      res_long_prop_trim <- residuals(lm(trimmed_data$x ~ -1 + trimmed_data$Z_conf + trimmed_data$Z_w))
      long_wghts_trim <- res_long_prop_trim / (crossprod(trimmed_data$x, res_long_prop_trim)[1, 1])
      long_wghts_trim_full <- ifelse(seq_len(n) %in% trimmed_data$trimmed_row_indices, 0, long_wghts_trim)
      rows_to_delete <- trimmed_data$trimmed_row_indices

      return(list(
        trim_reg = trim_reg,
        trimmed_data = trimmed_data,
        res_long_prop_trim = res_long_prop_trim,
        long_wghts_trim = long_wghts_trim,
        long_wghts_trim_full = long_wghts_trim_full,
        num_covar_trim = num_covar_trim,
        rows_to_delete = rows_to_delete,
        ridgereg = ridgereg
      ))
    }

    multi_results <- if (length(long_NAs) > 0) multi_comps() else list(
      ridgereg = NULL,
      trimmed_data = NULL,
      num_covar_trim = 0,
      rows_to_delete = integer(0)
    )
    
    shortreg <- lm(Y ~ -1 + x + Z_conf)
    num_covar_short <- ncol(Z_conf)

    residuals_df <- determine_residuals(long_NAs, long_reg, multi_results$ridgereg)
    sigma <- sqrt(sum(residuals_df$res^2) / residuals_df$df)
    wghts <- crossprod(short_wghts, Z_w)
    return(c(multi_results, covariate_comps,
            list(
               res_long_prop = res_long_prop, long_wghts = long_wghts,
               short_prop = short_prop, short_wghts = short_wghts, res_prop = res_prop, 
               resid_Z_w_Z_uc = resid_Z_w_Z_uc, Z_w_resid_Z_w_Z_uc = Z_w_resid_Z_w_Z_uc, Z_w_res_prop = Z_w_res_prop,
               B = B, inv_B = solve(B),
               Z_conf = Z_conf
             ),
             list(long_reg = long_reg, long_NAs = long_NAs, num_covar_long = num_covar_long, num_covar_trim = multi_results$num_covar_trim),
             list(shortreg = shortreg, num_covar_short = num_covar_short),
             list(residuals_df = residuals_df, sigma = sigma, wghts = wghts, rows_to_delete = multi_results$rows_to_delete)))
  })
}
