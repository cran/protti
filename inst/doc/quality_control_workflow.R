## ----include = FALSE----------------------------------------------------------
test_protti <- identical(Sys.getenv("TEST_PROTTI"), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, eval = test_protti, message = FALSE, warning = FALSE--------------
#  library(protti)
#  library(magrittr)
#  library(dplyr)

## ----create_synthetic_data, eval = test_protti--------------------------------
#  # by setting the seed we are making sure that the random object generation can be reproduced
#  set.seed(123)
#  
#  data <- create_synthetic_data(
#    n_proteins = 100,
#    frac_change = 0.05,
#    n_replicates = 3,
#    n_conditions = 2,
#    method = "effect_random",
#    additional_metadata = TRUE
#  )

## ----qc_cvs, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  input <- data %>%
#    # as the data is log2 transformed, we need to transform it back before calculating the CVs
#    mutate(raw_intensity = 2^peptide_intensity_missing)
#  
#  qc_cvs(
#    data = input,
#    grouping = peptide,
#    condition = condition,
#    intensity = raw_intensity,
#    plot = FALSE
#  )
#  
#  qc_cvs(
#    data = input,
#    grouping = peptide,
#    condition = condition,
#    intensity = raw_intensity,
#    plot = TRUE,
#    plot_style = "violin"
#  )

## ----qc_ids, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_ids(
#    data = input,
#    sample = sample,
#    grouping = protein,
#    intensity = peptide_intensity_missing,
#    condition = condition,
#    plot = FALSE
#  )
#  
#  qc_ids(
#    data = input,
#    sample = sample,
#    grouping = protein,
#    intensity = peptide_intensity_missing,
#    condition = condition,
#    title = "Protein identifications per sample",
#    plot = TRUE
#  )

## ----qc_peptide_type, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_peptide_type(
#    data = input,
#    sample = sample,
#    peptide = peptide,
#    pep_type = pep_type,
#    method = "intensity",
#    intensity = raw_intensity,
#    plot = TRUE,
#    interactive = FALSE
#  )
#  
#  qc_peptide_type(
#    data = input,
#    sample = sample,
#    peptide = peptide,
#    pep_type = pep_type,
#    method = "count",
#    plot = TRUE,
#    interactive = FALSE
#  )

## ----qc_intensity_distribution_boxplot, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_intensity_distribution(
#    data = input,
#    sample = sample,
#    grouping = peptide,
#    intensity_log2 = peptide_intensity_missing,
#    plot_style = "boxplot"
#  )
#  
#  qc_median_intensities(
#    data = input,
#    sample = sample,
#    grouping = peptide,
#    intensity = peptide_intensity_missing
#  )

## ----qc_charge_states, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_charge_states(
#    data = input,
#    sample = sample,
#    grouping = peptide,
#    charge_states = charge,
#    method = "intensity",
#    intensity = raw_intensity,
#    plot = TRUE
#  )

## ----qc_missed_cleavages, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_missed_cleavages(
#    data = input,
#    sample = sample,
#    grouping = peptide,
#    missed_cleavages = n_missed_cleavage,
#    method = "intensity",
#    intensity = raw_intensity,
#    plot = TRUE
#  )

## ----qc_sequence_coverage, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_sequence_coverage(
#    data = input,
#    protein_identifier = protein,
#    coverage = coverage
#  )

## ----qc_peak_width, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_peak_width(
#    data = input,
#    sample = sample,
#    intensity = peptide_intensity_missing,
#    retention_time = retention_time,
#    peak_width = peak_width
#  )

## ----qc_data_completeness, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_data_completeness(
#    data = input,
#    sample = sample,
#    grouping = peptide,
#    intensity = peptide_intensity_missing,
#    plot = TRUE
#  )

## ----qc_intensity_distribution_histogram, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_intensity_distribution(
#    data = input,
#    grouping = peptide,
#    intensity_log2 = peptide_intensity_missing,
#    plot_style = "histogram"
#  )

## ----qc_sample_correlation, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_sample_correlation(
#    data = input,
#    sample = sample,
#    grouping = peptide,
#    intensity_log2 = peptide_intensity_missing,
#    condition = condition,
#    interactive = FALSE
#  )

## ----qc_pca, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  qc_pca(
#    data = data,
#    sample = sample,
#    grouping = peptide,
#    intensity = peptide_intensity_missing,
#    condition = condition,
#    digestion = NULL,
#    plot_style = "scree"
#  )
#  
#  qc_pca(
#    data = data,
#    sample = sample,
#    grouping = peptide,
#    intensity = peptide_intensity_missing,
#    condition = condition,
#    components = c("PC1", "PC2"),
#    plot_style = "pca"
#  )

## ----qc_ranked_intensities, eval = test_protti, fig.width = 6, fig.height = 4, fig.align = "center"----
#  # Plot ranked peptide intensities
#  qc_ranked_intensities(
#    data = data,
#    sample = sample,
#    grouping = peptide,
#    intensity_log2 = peptide_intensity,
#    plot = TRUE,
#  )

