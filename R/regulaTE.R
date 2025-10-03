#' Function performing inference under a heterogeneity bound that end users call.
#'
#' @param outcome A numeric vector; realizations of the outcome variable.
#' @param treatment A numeric vector; realizations of the treatment variable.
#' @param covariates A numeric matrix; realizations of the control covariates.
#' @param confounders A numeric matrix; realizations of the confounders. Defaults to \code{NULL}.
#' If \code{NULL}, the function will use the covariates as confounders.
#' @param parameter A string; the parameter of interest, which must be \code{"ATE"}, \code{"ATT"}, or \code{"ATU"}. Defaults to \code{"ATE"}.
#' @param C A numeric scalar or vector; contains non-negative values specifying the heterogeneity bound or
#' \code{NA} if one wishes the program to estimate \code{C} using the data.
#' @param se A string; either \code{"het"}, or \code{"hom"}. \code{"hom"} computes homoskedastic standard errors. \code{"het"} computes 
#' heteroskedastic-robust standard errors. If \code{"het"} is chosen and \code{cluster} is specified, then it computes clustered
#'  heteroskedastic-robust standard errors.
#' @param cluster A vector; cluster identifier for computing clustered standard errors. Defaults to \code{NULL}. If \code{NULL}, non-clustered standard errors are used.
#' @param sig_level A numeric vector; significance level for confidence intervals. Defaults to 0.05.
#' @param df_corr indicating whether to apply degrees-of-freedom correction. Defaults to \code{TRUE}.
#' @param digits An integer vector; number of decimal places to use when formmatting numeric values. Defaults to 4.
#' @param trimmed_data A list; contains the matrices/vectors for the outcome (\code{Y}), 
#' treatment indicator (\code{x}), and interaction variables (\code{Z_w}), plus either 
#' covariates (\code{Z_uc}) or confounders (\code{Z_conf}). Used when trimming is necessary due to lack of 
#' overlap and the user wants to specify the trimmed data directly. Defaults to \code{NULL}, 
#' in which case the function trims the data automatically.
#' @param scaling A numeric scalar; a scaling factor for the maximum \eqn{\lambda} value. Defaults to 2.

#' @returns A list of two data frames:
#' \itemize{
#'   \item \code{summary}: A data frame; contains the sample size, number of covariates in the
#'   short and long regressions, and the Lindeberg weight.
#'   \item \code{estimates}: A data frame; contains point estimates, standard errors, and confidence
#'   intervals for each specification (Opt, Short, Short BC, and Long).
#' }
#'@export
#' 
regulaTE <- function(outcome, treatment, covariates, confounders = NULL, parameter = "ATE", C = NA,
                     se = "hom", cluster = NULL, sig_level = 0.95, df_corr = TRUE, digits = 4, trimmed_data = NULL, scaling = 2) {

  check_valid_inputs(outcome, treatment, covariates, parameter, C,
                                 se, cluster, sig_level, df_corr, digits, confounders, trimmed_data)

  covar_comps <- covariate_comps(treatment, covariates, parameter)

  reg_est <- reg_estimates(outcome, treatment, covar_comps, confounders, trimmed_data)

  if (is.null(cluster)&& se == "het") cluster <- seq_len(covar_comps$n)

  # Single C provided
  if (length(C) == 1) {
    output <- opt_het_bound(outcome, reg_est, cluster, C = C, alpha = 1 - sig_level, df_corr = df_corr, scaling = scaling)
    estimates_df <- clean_estimates_table(output$estimates_table, parameter, digits, sig_level)
    console_output(output, estimates_df, parameter, digits, sig_level, C, multi = FALSE)
    if (!is.na(C)) estimates_df$C <- C
    basic_info <- data.frame(
      sample_size = output$samp_size,
      num_covar_short = output$num_covar_short,
      num_covar_long = output$num_covar_long,
      lind_wt = output$lind_wt
    )
    return(invisible(list(summary = basic_info, estimates = estimates_df)))
  }

  # Multiple C's provided
  results <- lapply(C, function(current_C) {
    output <- opt_het_bound(outcome, reg_est, cluster, C = current_C, alpha = 1 - sig_level, df_corr = df_corr, scaling = scaling)
    estimates_df <- clean_estimates_table(output$estimates_table, parameter, digits, sig_level)

    estimates_df <- cbind(
      Specification = rownames(estimates_df),
      C = current_C,
      estimates_df,
      row.names = NULL
    )

    return(list(output = output, estimates_df = estimates_df))
  })

  stacked_estimates <- do.call(rbind, lapply(results, `[[`, "estimates_df"))
  rownames(stacked_estimates) <- NULL

  rep_output <- results[[if (any(sapply(C, is.na))) which(sapply(C, is.na))[1] else 1]]$output
  basic_info <- data.frame(
    sample_size = rep_output$samp_size,
    num_covar_short = rep_output$num_covar_short,
    num_covar_long = rep_output$num_covar_long
  )
  console_output(rep_output, NULL, parameter, digits, sig_level, C, multi = TRUE)
  return(invisible(list(summary = basic_info, estimates = stacked_estimates)))
}

#' Verify user inputs to regulaTE() are valid.
#'
#' @param outcome A numeric vector; realizations of the outcome variable.
#' @param treatment A numeric vector; realizations of the treatment variable.
#' @param covariates A numeric matrix; realizations of the control covariates.
#' @param confounders A numeric matrix; realizations of the confounders. Defaults to \code{NULL}.
#' If \code{NULL}, the function will use the covariates as confounders.
#' @param parameter A string; the parameter of interest, which must be \code{"ATE"}, \code{"ATT"}, or \code{"ATU"}. Defaults to \code{"ATE"}.
#' @param C A numeric scalar or vector; contains non-negative values specifying the heterogeneity bound or
#' \code{NA} if one wishes the program to estimate \code{C} using the data.
#' @param se A string; either \code{"het"}, or \code{"hom"}. \code{"hom"} computes homoskedastic standard errors. \code{"het"} computes 
#' heteroskedastic-robust standard errors. If \code{"het"} is chosen and \code{cluster} is specified, then it computes clustered
#'  heteroskedastic-robust standard errors.
#' @param cluster A vector; cluster identifier for computing clustered standard errors. Defaults to \code{NULL}. If \code{NULL}, non-clustered standard errors are used.
#' @param sig_level A numeric vector; significance level for confidence intervals. Defaults to 0.05.
#' @param df_corr indicating whether to apply degrees-of-freedom correction. Defaults to \code{TRUE}.
#' @param digits An integer vector; number of decimal places to use when formmatting numeric values. Defaults to 4.
#' @param trimmed_data A list; contains the matrices/vectors for the outcome (\code{Y}), 
#' treatment indicator (\code{x}), and interaction variables (\code{Z_w}), plus either 
#' covariates (\code{Z_uc}) or confounders (\code{Z_conf}). Used when trimming is necessary due to lack of 
#' overlap and the user wants to specify the trimmed data directly. Defaults to \code{NULL}, 
#' in which case the function trims the data automatically.
#' @returns None if inputs are valid; otherwise, the function stops the program.
#' 
check_valid_inputs <- function(outcome, treatment, covariates, parameter, C,
                               se, cluster, sig_level, df_corr, digits, confounders, trimmed_data) {
  errors <- c()

  if (!is.vector(outcome)) errors <- c(errors, "'outcome' must be a vector.")
  if (!is.vector(treatment)) errors <- c(errors, "'treatment' must be a vector.")
  if (!all(treatment %in% c(0, 1))) errors <- c(errors, "'treatment' must be a binary vector containing only 0s and 1s.")
  if (length(outcome) != length(treatment)) errors <- c(errors, "'outcome' and 'treatment' must be of the same length.")
  if (!is.matrix(covariates)) errors <- c(errors, "'covariates' must be a matrix.")
  if (nrow(covariates) != length(outcome)) errors <- c(errors, "The number of rows in 'covariates' must equal the length of 'outcome' and 'treatment'.")
  if (!parameter %in% c("ATE", "ATT", "ATU")) errors <- c(errors, "'parameter' must be 'ATE', 'ATT', or 'ATU'.")
  if (is.null(C)) errors <- c(errors, "No null 'C' values allowed.")
  if (!is.numeric(C) && !all(is.na(C))) errors <- c(errors, "'C' must be either NA or numeric.")
  if (is.numeric(C) && any(C < 0, na.rm = TRUE)) errors <- c(errors, "'C' must contain only non-negative values.")
  if (any(duplicated(C))) errors <- c(errors, "'C' must not contain duplicated values.")
  if (!is.numeric(sig_level) || length(sig_level) != 1) errors <- c(errors, "'sig_level' must be a single numeric value.")
  if (!is.logical(df_corr) || length(df_corr) != 1) errors <- c(errors, "'df_corr' must be TRUE or FALSE.")
  if (!is.numeric(digits) || length(digits) != 1 || digits <= 0 || digits != as.integer(digits)) errors <- c(errors, "'digits' must be a positive integer.")
  if (se == "hom" && !is.null(cluster)) errors <- c(errors, "For homoskedastic standard errors ('hom'), 'cluster' should not be specified.")
  if (se == "het" && !is.null(cluster) && length(cluster) != length(outcome)) errors <- c(errors, "'cluster' must be the same length as 'outcome' and 'treatment' if specified.")
  if (!is.null(confounders)){
    if (!is.matrix(confounders)) errors <- c(errors, "'confounders' must be a matrix if provided.")
    if (nrow(covariates) != length(outcome)) errors <- c(errors, "The number of rows in 'covariates' must equal the length of 'outcome' and 'treatment' if provided.")
  }
  if (!is.null(trimmed_data)) {
    if (!is.list(trimmed_data)) errors <- c(errors, "'trimmed_data' must be a list if provided, containing: Y (vector), x (vector), and Z_uc (matrix).")
    
    required_cols <- c("Y", "x", "Z_uc")
    missing_cols <- setdiff(required_cols, names(trimmed_data))
    if (length(missing_cols) > 0) {
      explanation <- paste0(
        "'trimmed_data' must contain the following named columns:\n",
        "- Y: the outcome variable\n",
        "- x: the treatment indicator\n",
        "- Z_uc: covariate matrix\n",
        "- Z_conf: confounder matrix (optional)\n",
      )
      errors <- c(errors, paste0("Missing columns in 'trimmed_data': ",
      paste(missing_cols, collapse = ","), "\n", explanation))
    }
    if (!is.null(cluster)){
      if (!"trimmed_row_indices" %in% names(trimmed_data)) errors <- c(errors, "For clustered SEs, 'trimmed_data$trimmed_row_indices' must be provided.")
      if (!(length(trimmed_data$trimmed_row_indices) + length(trimmed_data$Y) == length(outcome))) {
        errors <- c(errors, "The length of 'trimmed_data$trimmed_row_indices' does not match the size of the original data. These must satisfy: trimmed rows + remaining rows = total rows.")
      }
    }
  }


  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }
}
