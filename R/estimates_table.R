#' Quantile of folded standard normal distribution
#'
#' Computes the quantile of the folded standard normal distribution, i.e., the distribution of the absolute value
#' of a normally distributed variable with mean \eqn{x} and standard deviation \eqn{1}.
#' If \eqn{x > 3}, the function uses a normal approximation; otherwise, it uses the square root of a noncentral chi-squared quantile.
#'
#' @importFrom stats qnorm qchisq
#' 
#' @param x A numeric vector; the mean of the standard normal distribution.
#' @param q A numeric vector; the desired quantile level between \eqn{0} and \eqn{1}.
#'
#' @return A numeric vector; the \eqn{q}-th quantile of \eqn{|Z|}, where \eqn{Z \sim N(x, 1)}.
#' 
absnorm_q <- function(x, q) {
  if (x > 3) {
    cvalpha <- x + qnorm(q)
  } else {
    cvalpha <- sqrt(qchisq(q, df = 1, ncp = x ^ 2))
  }
  return(cvalpha)
}

#' Cluster-robust standard error computation
#'
#' Computes cluster-robust, heteroskedasticity-robust standard errors for a linear combination of residuals and weights,
#' with optional finite-sample degrees-of-freedom correction. It is thus equivalent to the \code{HC1} estimator (default in Stata).
#'
#' @param wghts A numeric vector; weights applied to the residuals.
#' @param residuals A numeric vector; residuals.
#' @param cluster A vector; cluster identifiers for each obsevation.
#' @param num_covar An integer vector; the number of covariates in the regression model.
#' @param df_corr A logical vector; if \code{TRUE}, applies a finite-sample correction to the standard error calculation.
#'
#' @return A numeric scalar corresponding to the cluster-robust standard error.
#'
cr_se_loop <- function(wghts, residuals, cluster, num_covar, df_corr = TRUE) {
  residuals <- as.vector(residuals)
  cluster <- as.factor(cluster)
  n <- length(residuals)
  G <- length(levels(cluster))

  if (G == n) {
    vhat <- sum((residuals * wghts)^2)
  } else {
    cluster_sums <- tapply(residuals * wghts, cluster, sum)
    vhat <- sum(cluster_sums^2)
  }

  if (df_corr && G > 1) {
    df_fac <- (G / (G - 1)) * ((n - 1) / (n - num_covar))
  } else {
    df_fac <- 1
  }
  return (sqrt(vhat * df_fac))
}
#' Generate a table of estimates for each regression specifcation
#'
#' @param opt_results A list; output from running \code{opt_lam()}.
#' @param reg_estimates A list; components obtained from running \code{reg_estimates()}.
#' @param sigma A numeric vector; residual standard deviation used to compute robust SEs when clustering is not specified.
#' @param alpha A numeric vector; between \code{0} and \code{1}; significance level for CIs.
#' @param zalpha A numeric vector; standard normal quantile for two-sided CI.
#' @param cva A numeric vector; quantile from folded normal distribution used for bias correction of the short regression.
#' @param max_bias_short A numeric vector; the maximum bias of the short regression.
#' @param max_bias_trim A numeric vector; the maximum bias of the trimmed long regression
#' @param cluster A vector; optional cluster identifiers for each obsevation. If \code{NULL}, standard errors are not clustered.
#' @param df_corr A logical vector; indicating whether to apply degrees-of-freedom correction. Defaults to \code{TRUE}.
#' 
#' @returns A numeric matrix; each row corresponds to a regression specification and includes
#' point estimates (PE), standard errors (SE), confidence intervals (CI LB and CI UB), and maximum bias, and critical value.
#'
generate_estimates_table <- function(opt_results, reg_estimates, sigma, alpha, zalpha, cva, max_bias_short, max_bias_trim, cluster, df_corr) {

  combine_estimates <- function(pe, se, cv, bias) {
    return (c(pe, se, pe - se*cv, pe+se*cv, bias, cv))
  }

  with(reg_estimates, {
    if (is.null(cluster)) {
      opt_estimates <- combine_estimates(opt_results$bhat, opt_results$se, opt_results$CIlen / opt_results$se, opt_results$max_bias)
      se_with_short_wghts <- sigma * sqrt(sum(short_wghts^2))
      if (length(long_NAs) > 0) {
        se_with_trim_wghts <- sigma * sqrt(sum(long_wghts_trim^2))
        long_estimates <- combine_estimates(trim_reg$coefficients[1], se_with_trim_wghts, zalpha, max_bias_trim)
        long_name <- "Trimmed Long"
      } else {
        se_with_long_wghts <- sigma * sqrt(sum(long_wghts^2))
        long_estimates <- combine_estimates(long_reg$coefficients[1], se_with_long_wghts, zalpha, 0)
        long_name <- "Long"
      }

      short_estimates <- combine_estimates(shortreg$coefficients[1], se_with_short_wghts, zalpha, max_bias_short)
      bc_short_estimates <- combine_estimates(shortreg$coefficients[1], se_with_short_wghts, cva, max_bias_short)

      estimates_table <- rbind(opt_estimates, short_estimates, bc_short_estimates, long_estimates)
      row.names(estimates_table) <- c("regulaTE", "Short", "Short BC", long_name)
    } else {
      if (length(long_NAs) > 0) {
        opt_se <- cr_se_loop(opt_results$wghts, residuals_df$res, cluster, num_covar_trim + 1, df_corr)
      } else {
        opt_se <- cr_se_loop(opt_results$wghts, residuals_df$res, cluster, num_covar_long + 1, df_corr)
      }
      opt_cva <- absnorm_q(opt_results$max_bias / opt_se, 1 - alpha)
      opt_estimates <- combine_estimates(opt_results$bhat, opt_se, opt_cva, opt_results$max_bias)

      short_se <- cr_se_loop(short_wghts, residuals_df$res, cluster, num_covar_long + 1, df_corr)
      short_estimates  <- combine_estimates(shortreg$coefficients[1], short_se, zalpha, max_bias_short)

      bc_cv_wghts  <- absnorm_q(max_bias_short / short_se,  1 - alpha)
      bc_estimates  <- combine_estimates(shortreg$coefficients[1], short_se, bc_cv_wghts, max_bias_short)

      if (length(long_NAs) > 0) {
        long_se <- cr_se_loop(long_wghts_trim, trim_reg$residuals, cluster[-rows_to_delete], num_covar_trim + 1, df_corr)
        long_estimates <- combine_estimates(trim_reg$coefficients[1], long_se, zalpha, max_bias_trim)
        long_name <- "Trimmed Long"
      } else {
        long_se <- cr_se_loop(long_wghts, residuals_df$res, cluster, num_covar_long + 1, df_corr)
        long_estimates <- combine_estimates(long_reg$coefficients[1], long_se, zalpha, 0)
        long_name <- "Long"
      }
      estimates_table <- rbind(opt_estimates, short_estimates, bc_estimates, long_estimates)
      row.names(estimates_table) <- c("regulaTE", "Short", "Short BC", long_name)
    }
    colnames(estimates_table) <- c("PE", "SE", "CI LB", "CI UB", "Max Bias", "Critical Value")
    return(estimates_table)
  })
}
