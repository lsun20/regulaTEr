#' Inference on treatment effects using a bound on heterogeneity
#'
#' Estimates the treatment effects (ATE, ATT, or ATU) using four different regressions: Opt, Short, Short BC, (Trimmed) Long.
#'
#' @param Y A numeric vector; outcome variable.
#' @param reg_estimates A list; output from \code{reg_estimates()}.
#' @param cluster A vector; cluster identifier for computing clustered standard errors. If \code{NULL}, non-clustered standard errors are used.
#' @param C A numeric vector; heterogeneity bound. If \code{NULL}, estimated using \code{cal_C()}.
#' @param alpha A numeric vector; significance level for confidence intervals. Defaults to 0.05.
#' @param df_corr A logical; indicating whether to apply degrees-of-freedom correction. Defaults to \code{TRUE}.
#' @param scaling A numeric scalar; a scaling factor for the maximum \eqn{\lambda} value. Defaults to 2.
#'
#' @returns A list with the following components:
#' \itemize{
#'   \item \code{estimates_table}: A numeric matrix; estimates for Opt, Short, Short BC, and (Trimmed) Long specs produced by \code{generate_estimates_table()}.
#'   \item \code{c_table}: A numeric matrix or \code{NA}; If \code{C} is estimated, a matrix containing PE, SE, and CI for \code{C} obtained from \code{generate_estimates_table()}; otherwise \code{NA}.
#'   \item \code{trimming_table}: A numeric matrix or \code{NA}; If trimming of the long regression is required, a matrix containing details of the trimming process obtained from \code{generate_estimates_table()}; otherwise \code{NA}.
#'   \item \code{lambda}: A numeric vector; value of the regularization parameter that minimizes CI length.
#'   \item \code{C_orig}:  A numeric vector or \code{NULL}; the user-supplied value of \code{C} or \code{NULL} if estimated using \code{cal_C()}.
#'   \item \code{samp_size}: A numeric vector; number of observations in the data.
#'   \item \code{num_covar_short}: A numeric vector; number of covariates used in the short regression.
#'   \item \code{num_covar_long}: A numeric vector; number of covariates used in the long regression.
#'   \item \code{lind_wt}: A numeric vector; the Lindeberg weight.
#'   \item \code{max_bias_short}: A numeric vector; the maximum bias of the short regression.
#' }
#'
opt_het_bound <- function(Y, reg_estimates, cluster, C, alpha = 0.05, df_corr = TRUE, scaling = 2) {
  with(reg_estimates, {
    C_orig <- C
    if (is.na(C)) {
      C_with_CI <- cal_C(Y, x, Z_uc, Z_w, B, trimmed_data, long_NAs)
      C <- sqrt(if (C_with_CI$C < 0) 0 else C_with_CI$C)
      c_table <- cbind(C_with_CI$C, C_with_CI$se, C_with_CI$lb, C_with_CI$ub)
      row.names(c_table) <- c(" ")
      colnames(c_table) <- c("C^2", "SE", "L", "U")
    } else {
      c_table <- NA
    }
    max_bias_short <- C * sqrt(crossprod(short_wghts, Z_w) %*% inv_B %*% t(crossprod(short_wghts, Z_w)))[1, 1]
    if (length(long_NAs) > 0) {
      max_bias_trim <- C * sqrt(crossprod(long_wghts_trim_full, Z_w) %*% inv_B %*% t(crossprod(long_wghts_trim_full, Z_w)))[1, 1]
    }

    cva <- absnorm_q(max_bias_short / (sqrt(sum(short_wghts ^ 2)) * sigma), 1 - alpha)
    half_length <- (sqrt(sum(short_wghts ^ 2)) * sigma) * cva

    lam_ran <- calculate_lambda_range(res_prop, B, resid_Z_w_Z_uc, x, short_wghts, scaling)

    opt_lam <- function(lambda, sigma, C) {
      delta <- delta <- solve(Z_w_resid_Z_w_Z_uc + lambda * B,  Z_w_res_prop)
      residuals <- res_prop - (resid_Z_w_Z_uc %*% delta)
      wghts <- residuals / sum(residuals * x)
      bhat <- sum(wghts * Y)
      pen <- sqrt(t(delta) %*% B %*% delta)[1, 1]
      max_bias <- C * crossprod(wghts, resid_Z_w_Z_uc) %*% delta / pen

      se <- sigma * sqrt(sum(wghts ^ 2))
      cvalpha <- absnorm_q(max_bias / se, 1 - alpha)

      CIlen <- se * cvalpha

      list(delta = delta, residuals = residuals, wghts = wghts, pen = pen, max_bias = max_bias,
           se = se, CIlen = CIlen, bhat = bhat, CI = c(bhat - CIlen, bhat + CIlen))
    }

    opt_obj <- optimise(f = function(lambda, ...) opt_lam(lambda, ...)$CIlen,
                        lam_ran, sigma = sigma, C = C)

    opt_results <- opt_lam(opt_obj$minimum, sigma, C)

    long_wghts_lind <- max((opt_results$wghts^2)/sum(opt_results$wghts^2))
    zalpha <- qnorm(1 - alpha / 2)

    estimates_table <- generate_estimates_table(opt_results, reg_estimates, sigma, alpha, zalpha, cva, max_bias_short, max_bias_trim, cluster, df_corr)
    return(list(
      estimates_table = estimates_table,
      c_table = c_table,
      trimming_table = trimmed_data$trimming_table,
      lambda = opt_obj$minimum,
      C_orig = C_orig,
      samp_size = n,
      num_covar_short = num_covar_short,
      num_covar_long = num_covar_long,
      lind_wt = long_wghts_lind,
      max_bias_short = max_bias_short
    ))
  })
}


#' Calculate an appropriate range for finding the optimal lambda value.
#' 
#' Determines a feasible search interval for the regularization parameter \eqn{\lambda} used in penalized propensity score regression.
#' The minimum value of this interval is always set to zero. 
#' The maximum is initialized using the largest \eqn{\lambda} value returned by \code{glmnet} because it recovers the fully regularized 
#' short regression coefficient estimates. The maximum \eqn{\lambda} is then adaptively increased by 2 until the \eqn{\ell_2}-distance 
#' between the weight vectors of this penalized propensity score regression and that of the fully penalized propensity score regression 
#' falls below a pre-specified tolerance.
#'
#' @importFrom stats lm predict
#'
#' @param res_prop A numeric vector; residuals from projecting \code{x} onto \code{Z_uc}.
#' @param B A numeric matrix; the covariance matrix of the demeaned \code{Z_uc}.
#' @param resid_Z_w_Z_uc A numeric matrix; residualized \code{Z_w} after projection onto \code{Z_uc}.
#' @param x A numeric vector; realizations corresponding to the treatment variable.
#' @param short_wghts A numeric vector; weights used in the short regression specification.
#' @param scaling A numeric scalar; a scaling factor for the maximum \eqn{\lambda} value. Defaults to 2.
#'
#' @returns A numeric vector; containing the minimum and maximum value of the lambda range.
#' 
calculate_lambda_range <- function(res_prop, B, resid_Z_w_Z_uc, x, short_wghts, scaling = 2) {
  Z_w_S <- resid_Z_w_Z_uc %*% solve(sqrt_mat(B))
  fit <- glmnet::glmnet(Z_w_S, res_prop, alpha = 0, standardize = FALSE)
  lambda_max <- max(fit$lambda)
  lambda_min <- 0

  get_l2_dist <- function(lambda) {
    fitted_vals <- predict(fit, newx = Z_w_S, s = lambda)
    residuals <- as.vector(res_prop - fitted_vals)
    wghts <- residuals / sum(residuals * x)
    return(sqrt(sum((short_wghts - wghts)^2)) / sqrt(sum(short_wghts^2)))
  }

  l2_dist <- get_l2_dist(lambda_max)
  while (l2_dist > 1e-8) {
    lambda_max <- lambda_max * 2
    l2_dist <- get_l2_dist(lambda_max)
  }
  lambda_max <- lambda_max * length(res_prop) * scaling

  return(c(lambda_min, lambda_max))
}
