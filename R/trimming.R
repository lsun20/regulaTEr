# Removes NA-causing rows (rows associated with limited overlap) in dataset.
#'
#' Helper function for \code{trim_data()}. Removes ("trims") all rows that have a value of 1 for any of the covariates
#' passed through in \code{na_covariates}.
#'
#' @param data A matrix; \code{cbind(Y, x, Z_uc)}, where \code{Y}, \code{x}, and \code{Z_uc} correspond to realizations of the outcome variable,
#' the treatment variable, and the covariate variables, respectively.
#' @param na_covariates A character vector; contains the names of interaction covariates with \code{NA} coefficients with the
#' \code{Z_w} prefix removed obtained from running \code{find_long_NAs()}.
#'
#' @returns A list containing the following components:
#' \itemize{
#'   \item \code{data}: The dataset with NA rows removed.
#'   \item \code{deleted_rows}: A numeric vector; contains the indices of the rows that were deleted.
#'   \item \code{rows_before}: A numeric scalar; equal to the number of rows in the data before trimming.
#'   \item \code{rows_after}: A numeric saclar; equal to the number of rows in the data after trimming.
#'   \item \code{rows_trimmed}: A numeric scalar; equal to the number of rows trimmed.
#'   \item \code{rows_trimmed_percentage}: A numeric scalar; equal to the percentage of how many rows were trimmed.
#' }
trim_na_rows <- function(data, na_covariates) {
  rows_before <- nrow(data)
  rowsToDelete <- unique(unlist(lapply(na_covariates, function(covar) which(data[, covar] == 1))))

  # cat("Breakdown of trimming by variable:\n")
  for (covar in na_covariates) {
    num_trimmed <- length(which(data[, covar] == 1))
    trimmed_percentage <- (num_trimmed / rows_before) * 100
    # cat(sprintf("%s: %d (%.2f%%)\n", covar, num_trimmed, trimmed_percentage))
  }

  data <- data[-rowsToDelete, ]
  rows_after <- nrow(data)
  rows_trimmed <- length(rowsToDelete)
  total_trimmed_percentage <- (rows_trimmed / rows_before) * 100

  # cat("Total rows before:", rows_before, "\n")
  # cat("Total rows after:", rows_after, "\n")
  # cat("Total rows trimmed:", rows_trimmed, "\n")
  # cat("Percent trimmed:", round(total_trimmed_percentage, 2), "%\n")

  list(data = data,
       deleted_rows = rowsToDelete,
       rows_before = rows_before,
       rows_after = rows_after,
       rows_trimmed = rows_trimmed,
       rows_trimmed_percentage = total_trimmed_percentage)
}

#' Removes NA-causing columns (columns associated limited overlap) in dataset.
#'
#' Helper function for \code{trim_data()}. Removes ("trims") all columns whose names are included in \code{na_covariates}.
#'
#' @param data A matrix; \code{cbind(Y, x, Z_uc)}, where \code{Y}, \code{x}, and \code{Z_uc} correspond to realizations of the outcome variable,
#' the treatment variable, and the covariate variables, respectively.
#' @param na_covariates A character vector; contains the names of interaction covariates with \code{NA} coefficients with the
#' \code{Z_w} prefix removed obtained from running \code{find_long_NAs()}.
#'
#' @returns A list containing the following components:
#' \itemize{
#'   \item \code{data}: The dataset with NA columns removed.
#'   \item \code{cols_before}: A numeric scalar; equal to the number of columns in the data before trimming.
#'   \item \code{cols_after}: A numeric saclar; equal to the number of columns in the data after trimming.
#'   \item \code{cols_trimmed}: A numeric scalar; equal to the number of columns trimmed.
#'   \item \code{cols_trimmed_percentage}: A numeric scalar; equal to the percentage of how many columns were trimmed.
#' }
trim_na_cols <- function(data, na_covariates) {
  cols_before <- length(colnames(data)) - 2 # subtract 2 to ignore outcome and treatment columns
  data <- data[, !colnames(data) %in% na_covariates, drop = FALSE]
  cols_after <- length(colnames(data)) - 2
  cols_trimmed <- cols_before - cols_after
  trimmed_percentage <- (cols_trimmed / cols_before) * 100

  # cat("Columns before:", cols_before, "\n")
  # cat("Columns after:", cols_after, "\n")
  # cat("Columns trimmed:", cols_trimmed, "\n")
  # cat("Percent trimmed:", round(trimmed_percentage, 2), "%\n")

  list(data = data,
       cols_before = cols_before,
       cols_after = cols_after,
       cols_trimmed = cols_trimmed,
       cols_trimmed_percentage = trimmed_percentage)
}


#' Trim data associated with limited overlap.
#'
#' Trims rows and columns that cause an undefined long regression. This function should be used only when \code{Z_conf = Z_uc} (i.e., no confounders).
#'
#' @param data A matrix; \code{cbind(Y, x, Z_uc)}, where \code{Y}, \code{x}, and \code{Z_uc} correspond to realizations of the outcome variable,
#' the treatment variable, and the covariate variables, respectively.
#' @param long_NAs A character vector; names of interaction covariates with \code{NA} coefficients, with the \code{Z_w} prefix removed. Possibly empty.
#' @param parameter A string; the parameter of interest, which must be \code{"ATE"}, \code{"ATT"}, or \code{"ATU"}. Defaults to \code{"ATE"}.
#'
#' @returns A list containing the following components:
#' \itemize{
#'   \item \code{Y}: A numeric vector; realizations of the outcome variable after trimming.
#'   \item \code{x}: A numeric vector; realizations of the treatment variable after trimming.
#'   \item \code{Z_uc}: A numeric matrix; realizations of the covariates after trimming.
#'   \item \code{Z_w}: A numeric matrix; realizations of the interacted covariates after trimming.
#'   \item \code{Z_tilde}: A numeric matrix; demeaned version of \code{Z_uc}, based on \code{parameter}.
#'   \item \code{B_trim}: A numeric matrix; covariance matrix of the demeaned \code{Z_uc}.
#'   \item \code{trimming_table}: A numeric matrix; 2-row summary of trimming across rows and columns.
#'   \item \code{trimmed_row_indices}: A numeric vector; indices of rows removed due to limited overlap.
#' }
#' Here "trimmed" means trimmed according to the procedures in \code{trim_na_rows()} and \code{trim_na_cols()}.
#' 
trim_data <- function(data, long_NAs, parameter) {
  # print("NAs in the long regression coefficients have been detected. Residual is from a regression with trimmed data.")

  trim_rows_result <- trim_na_rows(data, long_NAs)
  data <- trim_rows_result$data

  trim_cols_result <- trim_na_cols(data, long_NAs)
  data <- trim_cols_result$data

  y_trim <- data[, 1]
  x_trim <- data[, 2]
  Z_uc_trim <- data[, 3:ncol(data)]

  if (!all(Z_uc_trim[, ncol(Z_uc_trim)] == 1)) {
    Z_uc_trim <- cbind(Z_uc_trim, 1)
  }

  w <- switch(parameter,
              ATT = x_trim,
              ATU = 1 - x_trim,
              rep(1, length(x_trim)))  # default is ATE

  Z_tilde_trim <- apply(Z_uc_trim[, -ncol(Z_uc_trim)], 2, function(z) z - mean(w * z) / mean(w))
  Z_w_trim <- x_trim * Z_tilde_trim
  B_trim <- crossprod(Z_tilde_trim) / ncol(Z_uc_trim)

  trimming_table <- cbind(trim_rows_result$rows_before, trim_rows_result$rows_after,
                          trim_rows_result$rows_trimmed, trim_rows_result$rows_trimmed_percentage)
  trimming_table <- rbind(trimming_table, cbind(trim_cols_result$cols_before, trim_cols_result$cols_after,
                                                trim_cols_result$cols_trimmed, trim_cols_result$cols_trimmed_percentage))
  row.names(trimming_table) <- c("Trimming Breakdown for Rows", "Trimming Breakdown for Cols")
  colnames(trimming_table) <- c("Before", "After", "Trimmed", "Percent Trimmed")
  list(Y = y_trim,
       x = x_trim,
       Z_uc = Z_uc_trim,
       Z_w = Z_w_trim,
       Z_tilde = Z_tilde_trim,
       B_trim = B_trim,
       trimming_table = trimming_table,
       trimmed_row_indices = trim_rows_result$deleted_rows)
}
