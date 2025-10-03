#' Plots multiple estimates from regulaTE() with varying \code{C} values.
#'
#' @param estimates A data frame; \code{estimates_df} from \code{regulaTE()}.
#' @param y_axis_title A character string; the title for the y-axis. Defaults to "Outcome."
#' @param title A character string; the title for the plot. Defaults to "Estimates under Heterogeneity Bound for <Parameter>"
#' @param show_long A logical; whether to show the long regression results in the plot. Defaults to FALSE.
#' @returns Note. Only prints a plot.
#' 
#' @export
#' 
plot_regulaTE <- function(estimates, y_axis_title = "Outcome", title = NULL, show_long = FALSE) {
  check_valid_plot_inputs(estimates, y_axis_title, title)

  parameter <- names(estimates)[3]
  conf_int_col <- names(estimates)[6]

  if (is.null(title)) {
    title <- paste0("Estimates under Heterogeneity Bound for ", toupper(parameter))
  }

  df <- estimates
  df$PE <- df[[parameter]]

  conf_int <- strsplit(gsub("\\[|\\]", "", df[[conf_int_col]]), ",\\s*")
  df$L <- as.numeric(sapply(conf_int, `[`, 1))
  df$U <- as.numeric(sapply(conf_int, `[`, 2))
  legend_colors_base <- c(
    "Short CI" = "orange",
    "Short Point Estimate" = "orange",
    "Short BC Point Estimate" = "blue",
    "regulaTE Point Estimate" = "red",
    "Long CI" = "green",
    "Long Point Estimate" = "green"
  )
  legend_fills <- c(
    "regulaTE CI" = "lightpink",
    "Short BC CI" = "lightblue"
  )

  legend_order <- c(
    "regulaTE Point Estimate",
    "Short Point Estimate",
    "Short BC Point Estimate",
    "Long Point Estimate",
    "Short CI",
    "Short BC CI",
    "regulaTE CI",
    "Long CI")
  legend <- scale_color_manual(values = legend_colors_base)

  short_df <- df[df$Specification == "Short", ]
  if (nrow(short_df) > 0) {
    min_short_C <- min(short_df$C, na.rm = TRUE)
    short_point_estimate <- short_df$PE[short_df$C == min_short_C]
    short_ci_lower <- short_df$L[short_df$C == min_short_C]
    short_ci_upper <- short_df$U[short_df$C == min_short_C]
  } else {
    min_short_C <- NA
    short_point_estimate <- NA
    short_ci_lower <- NA
    short_ci_upper <- NA
  }

  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = min_short_C, y = short_ci_lower,
                                       xend = min_short_C, yend = short_ci_upper, color = "Short CI"),
                          linetype = "solid", size = 1, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(x = min_short_C, y = short_point_estimate, color = "Short Point Estimate"),
                        size = 3, na.rm = TRUE) +
    ggplot2::geom_ribbon(data = df[df$Specification == "Short BC", ],
                         ggplot2::aes(x = C, ymin = L, ymax = U, fill = "Short BC CI"),
                         alpha = 0.5, na.rm = TRUE) +
    ggplot2::geom_line(data = df[df$Specification == "Short BC", ],
                       ggplot2::aes(x = C, y = U), linetype = "solid", color = "blue",
                       show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_line(data = df[df$Specification == "Short BC", ],
                       ggplot2::aes(x = C, y = L), linetype = "solid", color = "blue",
                       show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_line(data = df[df$Specification == "Short BC", ],
                       ggplot2::aes(x = C, y = PE, color = "Short BC Point Estimate"),
                       linetype = "dashed", na.rm = TRUE) +
    ggplot2::geom_ribbon(data = df[df$Specification == "regulaTE", ],
                         ggplot2::aes(x = C, ymin = L, ymax = U, fill = "regulaTE CI"),
                         alpha = 0.5, na.rm = TRUE) +
    ggplot2::geom_line(data = df[df$Specification == "regulaTE", ],
                       ggplot2::aes(x = C, y = U), linetype = "solid", color = "red",
                       show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_line(data = df[df$Specification == "regulaTE", ],
                       ggplot2::aes(x = C, y = L), linetype = "solid", color = "red",
                       show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_line(data = df[df$Specification == "regulaTE", ],
                       ggplot2::aes(x = C, y = PE, color = "regulaTE Point Estimate"),
                       linetype = "dotted", na.rm = TRUE) +
    ggplot2::labs(
      title = title,
      x = "Heterogeneity Bound C",
      y = y_axis_title
    ) +
    ggplot2::scale_color_manual(values = legend_colors_base, breaks = legend_order) +
    ggplot2::scale_fill_manual(values = legend_fills, breaks = legend_order) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16, face = "bold")
    ) +
    scale_fill_manual(values = legend_fills) +
    legend +
    guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))



  if (show_long) {
    long_df <- df[df$Specification %in% c("Long","Trimmed Long"), ]
    max_long_C <- max(long_df$C, na.rm = TRUE)
    long_point_estimate <- long_df$PE[long_df$C == max_long_C]
    long_ci_lower <- long_df$L[long_df$C == max_long_C]
    long_ci_upper <- long_df$U[long_df$C == max_long_C]

    plot <- plot +
      ggplot2::geom_segment(ggplot2::aes(x = max_long_C, y = long_ci_lower,
                                         xend = max_long_C, yend = long_ci_upper, color = "Long CI"),
                            linetype = "solid", size = 1, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(x = max_long_C, y = long_point_estimate, color = "Long Point Estimate"),
                          size = 3, na.rm = TRUE)
  } 

  print(plot)
}

#' Checks if inputs into plot_regulaTE() are valid.
#'
#' @param estimates A data frame; \code{estimates_df} from regulaTE().
#' @param y_axis_title A character string; the title for the y-axis. Defaults to "Outcome."
#' @param title A character string; the tile for the graph. Defaults to "Estimates under Heterogeneity Bound for <Parameter>"
#'
#' @returns  None if inputs are valid; otherwise, the function stops the program.
#' 
check_valid_plot_inputs <- function(estimates, y_axis_title, title = NULL) {
  errors <- c()

  if (!is.data.frame(estimates)) {
    errors <- c(errors, "The 'estimates' argument must be a data frame.")
  } else {
    if (ncol(estimates) < 6) {
      errors <- c(errors, "The 'estimates' data frame must have at least 6 columns.")
    } else {
      if (names(estimates)[1] != "Specification") {
        errors <- c(errors, "The first column must be named 'Specification'.")
      }

      if (names(estimates)[2] != "C") {
        errors <- c(errors, "The second column must be named 'C'.")
      }

      parameter <- names(estimates)[3]
      if (!parameter %in% c("ATE", "ATU", "ATT")) {
        errors <- c(errors, "The third column must be one of 'ATE', 'ATU', or 'ATT'.")
      }

      if (names(estimates)[4] != "Std Err") {
        errors <- c(errors, "The fourth column must be named 'Std Err'.")
      }

      if (names(estimates)[5] != "Critical Value") {
        errors <- c(errors, "The fifth column must be named 'Critical Value'.")
      }

      conf_int_col <- names(estimates)[6]
      if (!grepl("^\\d+%\\s*Conf\\s*Int$", conf_int_col)) {
        errors <- c(errors, "The sixth column must be named in the format '<number>% Conf Int' (e.g., '95% Conf Int').")
      } else {
        conf_percentage <- as.numeric(gsub("%\\s*Conf\\s*Int", "", conf_int_col))
        if (is.na(conf_percentage) || conf_percentage < 0 || conf_percentage > 100) {
          errors <- c(errors, "The percentage in the sixth column name must be an integer between 0 and 100.")
        }
      }

      valid_specs <- c("regulaTE", "Short", "Short BC", "Long", "Trimmed Long")
      if (!all(estimates$Specification %in% valid_specs)) {
        invalid_specs <- unique(estimates$Specification[!estimates$Specification %in% valid_specs])
        errors <- c(errors, paste("The 'Specification' column contains invalid values:",
                                  paste(invalid_specs, collapse = ", "),
                                  ". Must be one of 'regulaTE', 'Short', 'Short BC', 'Long', or 'Trimmed Long'."))
      }

      for (col in 2:5) {
        if (!is.numeric(estimates[[col]])) {
          errors <- c(errors, paste("Column", col, " ('", names(estimates)[col], "') must be numeric."))
        }
      }

      conf_int_values <- as.character(estimates[[6]])
      invalid_format <- sapply(conf_int_values, function(ci) {
        if (!startsWith(ci, "[") || !endsWith(ci, "]") || (length(strsplit(ci, ",")[[1]]) != 2)) {
          return(TRUE)
        }
        no_brackets <- gsub("\\[|\\]", "", ci)
        parts <- strsplit(no_brackets, ",")[[1]]
        num1 <- suppressWarnings(as.numeric(parts[1]))
        num2 <- suppressWarnings(as.numeric(parts[2]))
        if (any(is.na(c(num1, num2))) || num1 > num2) {
          return(TRUE)
        }
        return(FALSE)
      })

      if (any(invalid_format)) {
        errors <- c(errors, paste0(
          "The confidence interval column (fifth column) must contain strings in the format '[lower,upper]', ",
          "where lower <= upper. Invalid entries found at rows: ", paste(which(invalid_format), collapse = ", ")
        ))
      }

    }

    if (!is.character(y_axis_title) || nchar(y_axis_title) == 0) {
      errors <- c(errors, "The 'y_axis_title' argument must be a non-empty string.")
    }

    if (!is.null(title) && (!is.character(title) || nchar(title) == 0)) {
      errors <- c(errors, "The 'title' argument must be NULL or a non-empty string.")
    }
  }

  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }
}
