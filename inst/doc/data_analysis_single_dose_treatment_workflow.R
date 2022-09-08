## ---- include = FALSE---------------------------------------------------------
test_protti <- identical(Sys.getenv("TEST_PROTTI"), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, eval = test_protti, message = FALSE, warning = FALSE--------------
library(protti)
library(magrittr)
library(dplyr)

## ----eval=FALSE---------------------------------------------------------------
#  # To read in your own data you can use read_protti()
#  your_data <- read_protti(filename = "mydata/data.csv")

## ----load_data, eval = test_protti--------------------------------------------
data("rapamycin_10uM")

## ----median_normalisation_filtering, eval = test_protti, message = FALSE, warning = FALSE----
data_normalised <- rapamycin_10uM %>%
  filter(eg_is_decoy == FALSE) %>%
  mutate(intensity_log2 = log2(fg_quantity)) %>%
  normalise(sample = r_file_name, 
            intensity_log2 = intensity_log2,
            method = "median")

data_filtered <- data_normalised %>%
  filter_cv(grouping = eg_precursor_id, 
            condition = r_condition, 
            log2_intensity = intensity_log2, 
            cv_limit = 0.25,
            min_conditions = 1)

## ----data_preparation_prot_pep, eval = test_protti, message = FALSE, warning = FALSE----
data_filtered_proteotypic <- data_filtered %>%
  filter(pep_is_proteotypic == TRUE)

## ----fetching_database_information, eval = test_protti, message = FALSE, warning = FALSE----
uniprot_ids <- unique(data_filtered_proteotypic$pg_protein_accessions)

uniprot <-
  fetch_uniprot(
    uniprot_ids = uniprot_ids,
    columns = c(
      "protein_name",
      "gene_names",
      "go_f",
      "xref_string",
      "cc_interaction",
      "ft_act_site",
      "ft_binding",
      "xref_pdb",
      "length",
      "sequence"
    )
  ) %>% 
  rename(pg_protein_accessions = accession)

data_filtered_uniprot <- data_filtered_proteotypic %>%
  left_join(y = uniprot,
            by = "pg_protein_accessions") %>%
  find_peptide(protein_sequence = sequence,
               peptide_sequence = pep_stripped_sequence) %>%
  assign_peptide_type(aa_before = aa_before,
               last_aa = last_aa, 
               aa_after = aa_after) %>%
  calculate_sequence_coverage(protein_sequence = sequence,
                    peptides = pep_stripped_sequence)

## ----coverage_plot, eval = test_protti, fig.align= "center", fig.width = 6, fig.height = 5----
qc_sequence_coverage(
  data = data_filtered_uniprot,
  protein_identifier = pg_protein_accessions,
  coverage = coverage
)

## ----t_test, eval = test_protti, message = FALSE, warning = FALSE-------------
diff_abundance_data <- data_filtered_uniprot %>%
  assign_missingness(
    sample = r_file_name,
    condition = r_condition,
    grouping = eg_precursor_id,
    intensity = normalised_intensity_log2,
    ref_condition = "control",
    completeness_MAR = 0.7,
    completeness_MNAR = 0.25,
    retain_columns = c(pg_protein_accessions, 
                       go_f, 
                       xref_string, 
                       start, 
                       end, 
                       length, 
                       coverage)
  ) %>%
  calculate_diff_abundance(
    sample = r_file_name,
    condition = r_condition,
    grouping = eg_precursor_id,
    intensity_log2 = normalised_intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    retain_columns = c(pg_protein_accessions, 
                       go_f, 
                       xref_string, 
                       start, 
                       end, 
                       length, 
                       coverage)
  ) 

## ----pval_distribution, eval = test_protti, message = FALSE, warning = FALSE, fig.align= "center", fig.width = 6, fig.height = 5----
pval_distribution_plot(data = diff_abundance_data,
                       grouping = eg_precursor_id,
                       pval = pval
                       )

## ----volcano_plot, eval = test_protti, fig.align= "center", fig.width = 6, fig.height = 5, message = FALSE, warning = FALSE----
volcano_plot(
  data = diff_abundance_data,
  grouping = eg_precursor_id,
  log2FC = diff,
  significance = pval,
  method = "target",
  target_column = pg_protein_accessions,
  target = "P62942",
  x_axis_label = "log2(fold change) Rapamycin treated vs. untreated",
  significance_cutoff = c(0.05, "adj_pval") 
)

# The significance_cutoff argument can also just be used for a 
# regular cutoff line by just providing the cutoff value, e.g.
# signficiance_cutoff = 0.05

## ----barcode_plot, eval = test_protti, fig.align = "center", fig.width = 6, message = FALSE, warning = FALSE----
FKBP12 <- diff_abundance_data %>%
  filter(pg_protein_accessions == "P62942")

barcode_plot(
  data = FKBP12,
  start_position = start,
  end_position = end,
  protein_length = length,
  coverage = coverage,
  colouring = diff,
  cutoffs = c(diff = 1, adj_pval = 0.05),
  protein_id = pg_protein_accessions
  )

## ----woods_plot, eval = test_protti, fig.align = "center", fig.width = 6, message = FALSE, warning = FALSE----

FKBP12 <- FKBP12 %>%
  mutate(significant = ifelse(adj_pval < 0.01, TRUE, FALSE))

woods_plot(
  data = FKBP12,
  fold_change = diff,
  start_position = start,
  end_position = end,
  protein_length = length,
  coverage = coverage,
  colouring = adj_pval,
  protein_id = pg_protein_accessions, 
  facet = FALSE,
  fold_change_cutoff = 1,
  highlight = significant
  )

## ----protile_plot, eval = test_protti, fig.align = "center", fig.width = 20, fig.height = 6, message = FALSE, warning = FALSE----
FKBP12_intensity <- data_filtered_uniprot %>% 
  filter(pg_protein_accessions == "P62942")

peptide_profile_plot(
  data = FKBP12_intensity,
  sample = r_file_name,
  peptide = eg_precursor_id,
  intensity_log2 = normalised_intensity_log2,
  grouping = pg_protein_accessions,
  targets = "P62942",
  protein_abundance_plot = FALSE
)

## ----additional_functions, eval=FALSE-----------------------------------------
#  diff_abundance_significant <- diff_abundance_data %>%
#    # mark significant peptides
#    mutate(is_significant = ifelse((adj_pval < 0.01 & abs(diff) > 1), TRUE, FALSE)) %>%
#    # mark true positive hits
#    mutate(binds_treatment = pg_protein_accessions == "P62942")
#  
#  ### GO enrichment using "molecular function" annotation from UniProt
#  
#  calculate_go_enrichment(
#    data = diff_abundance_significant,
#    protein_id = pg_protein_accessions,
#    is_significant = is_significant,
#    go_annotations_uniprot = go_f
#  )
#  
#  ### Network analysis
#  
#  network_input <- diff_abundance_significant %>%
#    filter(is_significant == TRUE)
#  
#  analyse_functional_network(data = network_input,
#                   protein_id = pg_protein_accessions,
#                   string_id = xref_string,
#                   binds_treatment = binds_treatment,
#                   organism_id = 9606)
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
#  diff_abundance_significant %>%
#    # columns containing proteins IDs are named differently
#    left_join(kegg, by = c("pg_protein_accessions" = "uniprot_id")) %>%
#    calculate_kegg_enrichment(protein_id = pg_protein_accessions,
#                    is_significant = is_significant,
#                    pathway_id = pathway_id,
#                    pathway_name = pathway_name)
#  
#  ### Treatment enrichment analysis
#  
#  calculate_treatment_enrichment(diff_abundance_significant,
#                       protein_id = pg_protein_accessions,
#                       is_significant = is_significant,
#                       binds_treatment = binds_treatment,
#                       treatment_name = "Rapamycin")
#  

