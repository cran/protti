#' Fitting four-parameter dose response curves (using parallel processing)
#'
#' This function is a wrapper around \code{fit_drc_4p} that allows the use of all system cores for
#' model fitting. It should only be used on systems that have enough memory available. Workers can
#' either be set up manually before running the function with \code{future::plan(multisession)} or
#' automatically by the function (maximum number of workers is 12 in this case). If workers are set
#' up manually the number of cores should be provided to \code{n_cores}. Worker can be terminated
#' after completion with \code{future::plan(sequential)}. It is not possible to export the
#' individual fit objects when using this function as compared to the non parallel function as
#' they are too large for efficient export from the workers.
#'
#' @details If data filtering options are selected, data is filtered based on multiple criteria.
#' In general, curves are only fitted if there are at least 5 conditions with data points present
#' to ensure that there is potential for a good curve fit. Therefore, this is also the case if no
#' filtering option is selected. Furthermore, a completeness cutoff is defined for filtering. By
#' default each entity (e.g. precursor) is filtered to contain at least 70% of total replicates
#' (adjusted downward) for at least 50% of all conditions (adjusted downward). This can be adjusted
#' with the according arguments. In addition to the completeness cutoff, also a significance cutoff
#' is applied. ANOVA is used to compute the statistical significance of the change for each entity.
#' The resulting p-value is adjusted using the Benjamini-Hochberg method and a cutoff of q <= 0.05
#' is applied. Curve fits that have a minimal value that is higher than the maximal value are
#' excluded as they were likely wrongly fitted. Curves with a correlation below 0.8 are not passing
#' the filtering. If a fit does not fulfill the significance or completeness cutoff, it has a
#' chance to still be considered if half of its values (+/-1 value) pass the replicate completeness
#' criteria and half do not pass it. In order to fall into this category, the values that fulfill
#' the completeness cutoff and the ones that do not fulfill it need to be consecutive, meaning
#' located next to each other based on their concentration values. Furthermore, the values that
#' do not pass the completeness cutoff need to be lower in intensity. Lastly, the difference
#' between the two groups is tested for statistical significance using a Welch's t-test and a
#' cutoff of p <= 0.1 (we want to mainly discard curves that falsly fit the other criteria
#' but that have clearly non-significant differences in mean). This allows curves to be considered
#' that have missing values in half of their observations due to a decrease in intensity. It can
#' be thought of as conditions that are missing not at random (MNAR). It is often the case that
#' those entities do not have a significant p-value since half of their conditions are not
#' considered due to data missingness.
#'
#' The final filtered list is ranked based on a score calculated on entities that pass the filter.
#' The score is the negative log10 of the adjusted ANOVA p-value scaled between 0 and 1 and the
#' correlation scaled between 0 and 1 summed up and divided by 2. Thus, the highest score an
#' entity can have is 1 with both the highest correlation and adjusted p-value. The rank is
#' corresponding to this score. Please note, that entities with MNAR conditions might have a
#' lower score due to the missing or non-significant ANOVA p-value. You should have a look at
#' curves that are TRUE for \code{dose_MNAR} in more detail.
#'
#' @param data a data frame that contains at least the input variables.
#' @param sample  a character column in the \code{data} data frame that contains the sample names.
#' @param grouping  a character column in the \code{data} data frame that contains the precursor,
#' peptide or protein identifiers.
#' @param response  a numeric column in the \code{data} data frame that contains the response
#' values, e.g. log2 transformed intensities.
#' @param dose  a numeric column in the \code{data} data frame that contains the dose values, e.g.
#' the treatment concentrations.
#' @param filter a character value that determines if models should be filtered and if they should
#' be filtered. The option \code{"pre"} is not available for parallel fitting of models. This is
#' because ANOVA adjusted p-values would be wrongly calculated because the dataset is split onto
#' multiple cores. Default is "post" and we recommend always using "post" because compared to
#' "none" only some additional columns are added that contain the filter information. For ANOVA
#' an adjusted p-value of 0.05 is used as a cutoff.
#' @param replicate_completeness a numeric value which similar to \code{completenss_MAR} of the
#' \code{assign_missingness} function sets a threshold for the completeness of data. In contrast
#' to \code{assign_missingness} it only determines the completeness for one condition and not the
#' comparison of two conditions. The threshold is used to calculate a minimal degree of data
#' completeness. The value provided to this argument has to be between 0 and 1, default is 0.7.
#' It is multiplied with the number of replicates and then adjusted downward. The resulting number
#' is the minimal number of observations that a condition needs to have to be considered "complete
#' enough" for the \code{condition_completeness} argument.
#' @param condition_completeness a numeric value which determines how many conditions need to at
#' least fulfill the "complete enough" criteria set with \code{replicate_completeness}. The
#' value provided to this argument has to be between 0 and 1, default is 0.5. It is multiplied with
#' the number of conditions and then adjusted downward. The resulting number is the minimal number
#' of conditions that need to fulfill the \code{replicate_completeness} argument for a peptide to
#' pass the filtering.
#' @param correlation_cutoff a numeric vector that specifies the correlation cutoff used for data
#' filtering.
#' @param log_logarithmic a logical value that indicates if a logarithmic or log-logarithmic model
#' is fitted. If response values form a symmetric curve for non-log transformed dose values, a
#' logarithmic model instead of a log-logarithmic model should be used. Usually biological dose
#' response data has a log-logarithmic distribution, which is the reason this is the default.
#' Log-logarithmic models are symmetric if dose values are log transformed.
#' @param retain_columns a vector that specifies columns that should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#' @param n_cores optional, a numeric value that specifies the number of cores used if workers
#' are set up manually.
#'
#' @return A data frame is returned that contains correlations of predicted to measured values as
#' a measure of the goodness of the curve fit, an associated p-value and the four parameters of
#' the model for each group. Furthermore, input data for plots is returned in the columns \code{plot_curve}
#' (curve and confidence interval) and \code{plot_points} (measured points).
#'
#' @import dplyr
#' @importFrom stats p.adjust
#' @importFrom rlang .data as_name enquo
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Load libraries
#' library(dplyr)
#'
#' set.seed(123) # Makes example reproducible
#'
#' # Create example data
#' data <- create_synthetic_data(
#'   n_proteins = 2,
#'   frac_change = 1,
#'   n_replicates = 3,
#'   n_conditions = 8,
#'   method = "dose_response",
#'   concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
#'   additional_metadata = FALSE
#' )
#'
#' # Perform dose response curve fit
#' drc_fit <- parallel_fit_drc_4p(
#'   data = data,
#'   sample = sample,
#'   grouping = peptide,
#'   response = peptide_intensity_missing,
#'   dose = concentration,
#'   retain_columns = c(protein, change_peptide)
#' )
#'
#' glimpse(drc_fit)
#'
#' head(drc_fit, n = 10)
#' }
parallel_fit_drc_4p <- function(data,
                                sample,
                                grouping,
                                response,
                                dose,
                                filter = "post",
                                replicate_completeness = 0.7,
                                condition_completeness = 0.5,
                                correlation_cutoff = 0.8,
                                log_logarithmic = TRUE,
                                retain_columns = NULL,
                                n_cores = NULL) {
  dependency_test <- c(
    furrr = !requireNamespace("furrr", quietly = TRUE),
    future = !requireNamespace("future", quietly = TRUE),
    parallel = !requireNamespace("parallel", quietly = TRUE),
    drc = !requireNamespace("drc", quietly = TRUE)
  )
  if (any(dependency_test)) {
    dependency_name <- names(dependency_test[dependency_test == TRUE])
    if (length(dependency_name) == 1) {
      message("Package \"",
        paste(dependency_name),
        "\" is needed for this function to work. Please install it.",
        call. = FALSE
      )
      return(invisible(NULL))
    } else {
      message("Packages \"",
        paste(dependency_name, collapse = "\" and \""),
        "\" are needed for this function to work. Please install them.",
        call. = FALSE
      )
      return(invisible(NULL))
    }
  }
  if (filter == "pre") {
    stop(strwrap('"pre" cannot be selected as a filter option for parallel
fitting. Use "post" or use the fit_drc_4p function.',
      prefix = "\n", initial = ""
    ))
  }
  . <- NULL
  terminate <- FALSE

  if (missing(n_cores)) {
    message("Setting up workers ... ", appendLF = FALSE)
    terminate <- TRUE
    n_cores <- parallel::detectCores()
    n_cores <- ifelse(n_cores > 12, 12, n_cores)
    oplan <- future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(oplan), add = TRUE)
    message("DONE", appendLF = TRUE)
  }

  pieces <- rep(
    1:n_cores,
    round(
      length(
        unique(dplyr::pull(data, {{ grouping }}))
      ) / n_cores
    ) + 1
  )[1:length(unique(dplyr::pull(data, {{ grouping }})))]

  pieces_mapping <- tibble::tibble({{ grouping }} := unique(dplyr::pull(data, {{ grouping }})), piece = pieces)

  input <- data %>%
    dplyr::left_join(pieces_mapping, by = rlang::as_name(rlang::enquo(grouping))) %>%
    split(dplyr::pull(., .data$piece))

  message("Performing model fit (this may take a while) ... ", appendLF = FALSE)

  result <- furrr::future_map_dfr(
    .x = input,
    .f = ~ protti::fit_drc_4p(.x,
      sample = {{ sample }},
      grouping = {{ grouping }},
      response = {{ response }},
      dose = {{ dose }},
      filter = filter,
      replicate_completeness = replicate_completeness,
      condition_completeness = condition_completeness,
      correlation_cutoff = correlation_cutoff,
      log_logarithmic = log_logarithmic,
      retain_columns = {{ retain_columns }},
      include_models = FALSE
    ),
    .options = furrr::furrr_options(globals = FALSE)
  )

  message("DONE", appendLF = TRUE)

  if (terminate == TRUE) {
    future::plan(future::sequential)
  }

  result <- result %>%
    dplyr::arrange(desc(.data$correlation))

  if (filter == "post") {
    result <- result %>%
      dplyr::mutate(anova_adj_pval = stats::p.adjust(.data$anova_pval, method = "BH")) %>%
      dplyr::mutate(anova_significant = ifelse(.data$anova_adj_pval > 0.05 | is.na(.data$anova_adj_pval),
        FALSE,
        TRUE
      )) %>%
      dplyr::mutate(passed_filter = ((.data$enough_conditions == TRUE &
        .data$anova_significant == TRUE) |
        .data$dose_MNAR == TRUE) &
        .data$correlation >= correlation_cutoff &
        .data$min_model < .data$max_model) %>%
      dplyr::group_by(.data$passed_filter) %>%
      dplyr::mutate(score = ifelse(.data$passed_filter,
        (scale_protti(-log10(.data$anova_pval), method = "01") + scale_protti(.data$correlation, method = "01")) / 2,
        NA
      )) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(.data$correlation)) %>%
      dplyr::arrange(dplyr::desc(.data$score)) %>%
      dplyr::select(-.data$rank) %>%
      tibble::rownames_to_column(var = "rank") %>%
      dplyr::mutate(rank = ifelse(!is.na(.data$score), as.numeric(.data$rank), NA))
  }

  return(result)
}
