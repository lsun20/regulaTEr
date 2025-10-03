#' Console output for regulaTE().
#'
#' Prints formatted console output for results from \code{opt_het_bound()}.
#' 
#' @param output A list; output from \code{opt_het_bound()}
#' @param estimates_table A data frame; produced by \code{produce_estimates_table()}.
#' @param parameter A character string; the parameter of interest, which must be one of \code{"ATE"}, \code{"ATT"}, or \code{"ATU"}.
#' @param digits An numeric vector; number of decimal places to display in the output.
#' @param sig_level A numeric vector; significance level for confidence intervals.
#' @param C A numeric vector; non-negative values specifying the heterogeneity bound or \code{NA} if one wishes the program to estimate \code{C} using the data.
#' @param multi A logical vector; whether the output is based on multiple \code{C} values. Defaults to \code{FALSE}.
#'
#' @returns None. Only prints to console.
console_output <- function(output, estimates_table, parameter, digits, sig_level, C, multi = FALSE) {
  cat("Parameter estimated:", parameter, "\n")
  cat("Sample size:", output$samp_size, "\n")
  cat("No. of covariates in short regression:", output$num_covar_short, "\n")
  cat("No. of covariates in long regression:", output$num_covar_long, "\n")

  if(!multi) {
    if (!is.na(C)) {
      cat("C (given by user):", C, "\n")
    } else {
      cat("C squared (estimated from data):", output$c_table[1, 1], "\n")
    }
    cat("Lindeberg weights:", output$lind_wt, "\n")
    cat("\n")
    cat("Estimates:\n")
    print_pretty_table(estimates_table, digits = digits)
  } else {
    cat("C values used:", paste(C, collapse = ", "), "\n")
  }
  if (!is.null(output$trimming_table)) {
    cat("\n")
    cat("Due to lack of overlap, the dataset has been trimmed accordingly:\n")
    cat("\n")
    trimming_table <- as.data.frame(output$trimming_table)
    trimming_table[] <- lapply(trimming_table, function(col) if (is.numeric(col)) round(col, digits) else col)
    rownames(trimming_table) <- c("Observations", "Covariates")
    print_pretty_table(trimming_table, digits = digits)
    cat("\n")
    cat("Note that one must therefore fix a seed to ensure replicability.")
  }
  if (any(is.na(C))) {
    cat("\n")
    c_table <- produce_c_table(as.data.frame(output$c_table), digits, sig_level)
    print_pretty_table(c_table, digits = digits)
  }
}

#' Prints a well-formatted table to the console 
#' 
#' Prints a well-formatted table to the console for a data frame of parameter or \code{C} estimates.
#'
#' @param df A data frame; contains parameter of \code{C} estimates from \code{produce_estimates_table()} or \code{produce_c_table()}.
#' @param digits An integer vector; number of decimal places to use when formmatting numeric values.
#' @param pad An integer vector; number of spaces between the columns. Defaults to 2.
#'
#' @returns None. Only prints to console.
#' 
print_pretty_table <- function(df, digits, pad = 2) {
  if (!is.null(rownames(df)) && any(rownames(df) != "")) {
    rownames_col <- format(rownames(df), justify = "left")
    df <- cbind(" " = rownames_col, df)
  }

  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      sapply(col, function(x) {
        if (is.na(x)) return("NA")
        format(x, justify = "right")
      })
    } else {
      format(col, justify = "centre")
    }
  })

  col_widths <- mapply(function(col, name) {
    max(nchar(as.character(col)), nchar(name))
  }, df, names(df))

  centered_headers <- mapply(function(name, width) {
    pad_total <- max(0, width - nchar(name))
    left_pad <- floor(pad_total / 2)
    right_pad <- ceiling(pad_total / 2)
    paste0(strrep(" ", left_pad), name, strrep(" ", right_pad))
  }, names(df), col_widths, SIMPLIFY = TRUE)

  cat(paste(centered_headers, collapse = strrep(" ", pad)), "\n")
  cat(paste(mapply(strrep, "-", col_widths), collapse = strrep(" ", pad)), "\n")

  for (i in 1:nrow(df)) {
    row <- paste(mapply(function(x, w) format(x, width = w, justify = "right"), df[i, ], col_widths), collapse = strrep(" ", pad))
    cat(row, "\n")
  }
}


#' Clean parameter estimates table for formatted output.
#' 
#' Produces a cleaned version of the estimates table from \code{opt_het_bound()}, formatted for display via \code{print_pretty_table()}.
#'
#' @param table A data frame; contains the point estimates, standard errors, and confidence intervals corresponding to the parameter estimates
#' from \code{opt_het_bound()}.
#' @param parameter A character vector; the parameter of interest, which must be one of \code{"ATE"}, \code{"ATT"}, or \code{"ATU"}.
#' @param digits An integer vector; number of decimal places to use when formmatting numeric values.
#' @param sig_level A numeric vector; significance level for confidence intervals.
#'
#' @returns A data frame; contains the point estimates, standard errors, and confidence intervals corresponding
#' to the parameter estimates from \code{opt_het_bound()} formalized appropriately for console output.
clean_estimates_table <- function(table, parameter, digits, sig_level) {
  ci_strings <- paste0("[", format(round(pmin(table[, 3], table[, 4]), digits), nsmall = digits), ", ",
                       format(round(pmax(table[, 3], table[, 4]), digits), nsmall = digits), "]")

  ci_colname <- paste0(round(100 * sig_level), "% Conf Int")

  final_table <- data.frame(
    ` ` = round(table[, 1, drop = FALSE], digits),
    `Std Err` = round(table[, 2], digits),
    `Max Bias` = round(table[, 5], digits),
    `Critical Value` = round(table[, 6], digits),
    check.names = FALSE,
    row.names = rownames(table)
  )

  colnames(final_table)[1] <- parameter
  final_table[[ci_colname]] <- ci_strings
  final_table <- final_table[, c(1, 2, 4, 5, 3)]

  return(final_table)
}

#' Produces \code{C} estimates table for formatted output.
#' 
#' Produces a cleaned data frame containing the estimated \eqn{C^2}, standard error, and confidence interval, formatted for display via \code{print_pretty_table()}.
#'
#' @param table A data frame; contains the estimates.
#' @param digits An integer vector; number of decimal places to use when formmatting numeric values.
#' @param sig_level A numeric vector; significance level for confidence intervals.
#'
#' @returns A data frame; contains the point estimate, standard error, and confidence interval corresponding to the
#' estimated \eqn{\vec{C}} from \code{opt_het_bound()} formalized appropriately for console output.
produce_c_table <- function(table, digits, sig_level) {
  cat("\n")
  cat("Since C contains an NA value, the program has estimated C:\n")
  cat("\n")
  row_name <- rownames(table)[1]

  ci_string <- paste0("[", format(round(min(table[, 3], table[, 4]), digits), nsmall = digits),
                      ", ", format(round(max(table[, 3], table[, 4]), digits), nsmall = digits), "]")

  final_table <- data.frame(
    `C^2` = round(table[, 1], digits),
    `Std Err` = round(table[, 2], digits),
    check.names = FALSE,
    row.names = row_name
  )
  final_table[[paste0(round(100 * sig_level), "% Conf Int")]] <- ci_string
  return(final_table)
}
