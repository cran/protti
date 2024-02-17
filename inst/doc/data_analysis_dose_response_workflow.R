## ----include = FALSE----------------------------------------------------------
build_vignette_on_cran <- identical(Sys.getenv("BUILD_VIGNETTE"), "true")
test_protti <- identical(Sys.getenv("TEST_PROTTI"), "true") & build_vignette_on_cran
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----CRAN_comment, message=FALSE, warning=FALSE, echo=FALSE-------------------
if (build_vignette_on_cran == FALSE) {
  print("!!! IMPORTANT !!!")
  print("This Vignette has not been built completely on CRAN due to size limitations.")
  print("Please check the correct version here: ")
  print("https://jpquast.github.io/protti/articles/https://jpquast.github.io/protti/articles/data_analysis_dose_response_workflow.html")
}

## ----load_packages, eval = test_protti, warning = FALSE, message = FALSE------
#  # Load packages
#  library(protti)
#  library(dplyr)
#  library(magrittr)
#  library(ggplot2)

## ----use_read_protti, eval=FALSE----------------------------------------------
#  # Load data
#  rapamycin_dose_response <- read_protti("your_data.csv")

## ----load_data, eval = test_protti, warning = FALSE---------------------------
#  utils::data("rapamycin_dose_response")

## ----filter_transform_normalise, eval = test_protti---------------------------
#  # Filter, log2 transform and normalise data
#  data_normalised <- rapamycin_dose_response %>%
#    filter(eg_is_decoy == FALSE) %>%
#    mutate(intensity_log2 = log2(fg_quantity)) %>%
#    normalise(
#      sample = r_file_name,
#      intensity_log2 = intensity_log2,
#      method = "median"
#    ) %>%
#    filter(pep_is_proteotypic == TRUE)

## ----intensity_distribution, eval = test_protti, fig.align = "center", fig.width = 7, fig.height = 5----
#  qc_intensity_distribution(
#    data = data_normalised,
#    grouping = eg_precursor_id,
#    intensity = normalised_intensity_log2,
#    plot_style = "histogram"
#  )

## ----intensity_filtering, eval = test_protti----------------------------------
#  data_normalised <- data_normalised %>%
#    filter(normalised_intensity_log2 > 5)

## ----sample_correlation, eval = test_protti, fig.align = "center", fig.width = 7, fig.height = 5----
#  qc_sample_correlation(
#    data = data_normalised,
#    sample = r_file_name,
#    grouping = eg_precursor_id,
#    intensity_log2 = normalised_intensity_log2,
#    condition = r_condition
#  )

## ----sample_pca, eval = test_protti, fig.align = "center", fig.width = 7, fig.height = 5, warning = FALSE, message = FALSE----
#  qc_pca(
#    data = data_normalised,
#    sample = r_file_name,
#    grouping = eg_precursor_id,
#    intensity = normalised_intensity_log2,
#    condition = r_condition
#  )

## ----model_fit, eval = test_protti--------------------------------------------
#  fit <- data_normalised %>%
#    fit_drc_4p(
#      sample = r_file_name,
#      grouping = eg_precursor_id,
#      response = normalised_intensity_log2,
#      dose = r_condition,
#      filter = "post",
#      retain_columns = c(pg_protein_accessions)
#    )
#  
#  # make sure to retain columns that you need later but that are not part of the function

## ----parallel_model_fit, eval = FALSE-----------------------------------------
#  # setup of cores. Make sure you have the future package installed
#  future::plan(future::multisession, workers = 3)
#  
#  # fit models in parallel
#  parallel_fit <- data_normalised %>%
#    parallel_fit_drc_4p(
#      sample = r_file_name,
#      grouping = eg_precursor_id,
#      response = normalised_intensity_log2,
#      dose = r_condition,
#      retain_columns = c(pg_protein_accessions),
#      n_cores = 3
#    )
#  
#  # remove workers again after you are done
#  future::plan(future::sequential)

## ----result_analysis, eval = test_protti, echo=FALSE, results='asis'----------
#  fit %>%
#    filter(rank <= 20) %>%
#    select(
#      rank,
#      score,
#      eg_precursor_id,
#      pg_protein_accessions,
#      anova_adj_pval,
#      correlation,
#      ec_50
#    ) %>%
#    mutate(
#      anova_adj_pval = format(anova_adj_pval, digits = 3),
#      correlation = format(correlation, digits = 3),
#      ec_50 = format(ec_50, digits = 2),
#      score = format(score, digits = 3)
#    ) %>%
#    knitr::kable(caption = "All hits")

## ----model_plot, eval = test_protti, fig.align = "center", fig.width = 7, fig.height = 5, message = FALSE, warning = FALSE----
#  # Model plotting
#  drc_4p_plot(fit,
#    grouping = eg_precursor_id,
#    dose = r_condition,
#    response = normalised_intensity_log2,
#    targets = "_VFDVELLKLE_.2",
#    unit = "pM",
#    export = FALSE
#  )

## ----uniprot_annotation, eval = test_protti-----------------------------------
#  # fetching of UniProt information
#  unis <- unique(fit$pg_protein_accessions)
#  uniprot <- fetch_uniprot(unis)
#  
#  # annotation of fit data based on information from UniProt
#  fit_annotated <- fit %>%
#    # columns containing proteins IDs are named differently
#    left_join(uniprot, by = c("pg_protein_accessions" = "accession")) %>%
#    # mark peptides that pass the filtering
#    mutate(passed_filter = !is.na(rank)) %>%
#    # create new column with prior knowledge about binding partners of treatment
#    mutate(binds_treatment = pg_protein_accessions == "P62942")

## ----post_analysis, eval = FALSE----------------------------------------------
#  ### GO enrichment using "molecular function" annotation from UniProt
#  
#  calculate_go_enrichment(fit_annotated,
#    protein_id = pg_protein_accessions,
#    is_significant = passed_filter,
#    go_annotations_uniprot = go_f
#  ) # column obtained from UniProt
#  
#  ### KEGG pathway enrichment
#  
#  # First you need to load KEGG pathway annotations from the KEGG database
#  # for your specific organism of interest. In this case HeLa cells were
#  # used, therefore the organism of interest is homo sapiens (hsa)
#  
#  kegg <- fetch_kegg(species = "hsa")
#  
#  # Next we need to annotate our data with KEGG pathway IDs and perform enrichment analysis
#  
#  fit %>%
#    # columns containing proteins IDs are named differently
#    left_join(kegg, by = c("pg_protein_accessions" = "uniprot_id")) %>%
#    calculate_kegg_enrichment(
#      protein_id = pg_protein_accessions,
#      is_significant = passed_filter,
#      pathway_id = pathway_id, # column name from kegg data frame
#      pathway_name = pathway_name
#    ) # column name from kegg data frame
#  
#  ### Treatment enrichment analysis
#  
#  calculate_treatment_enrichment(fit_annotated,
#    protein_id = pg_protein_accessions,
#    is_significant = passed_filter,
#    binds_treatment = binds_treatment,
#    treatment_name = "Rapamycin"
#  )
#  
#  ### Network analysis
#  
#  fit_annotated %>%
#    filter(passed_filter == TRUE) %>% # only analyse hits that were significant
#    analyse_functional_network(
#      protein_id = pg_protein_accessions,
#      string_id = xref_string, # column from UniProt containing STRING IDs
#      organism_id = 9606,
#      # tax ID can be found in function documentation or STRING database
#      binds_treatment = binds_treatment,
#      plot = TRUE
#    )

