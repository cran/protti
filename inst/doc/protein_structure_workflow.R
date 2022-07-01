## ---- include = FALSE---------------------------------------------------------
build_vignette_on_cran <- identical(Sys.getenv("BUILD_VIGNETTE"), "true")
test_protti <- identical(Sys.getenv("TEST_PROTTI"), "true") & build_vignette_on_cran
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)

## ----CRAN_comment, message=FALSE, warning=FALSE, echo=FALSE-------------------
if (build_vignette_on_cran == FALSE){
  print("!!! IMPORTANT !!!")
  print("This Vignette has not been built completely on CRAN due to size limitations.")
  print("Please check the correct version here: ")
  print("https://jpquast.github.io/protti/articles/protein_structure_workflow.html")
}

## ----load_packages, eval = test_protti, warning = FALSE, message = FALSE------
#  # Load packages
#  library(protti)
#  library(dplyr)
#  library(magrittr)
#  library(stringr)
#  library(tidyr)
#  library(ggplot2)

## ----load_data, eval = test_protti, warning = FALSE---------------------------
#  utils::data("ptsi_pgk")

## ----use_read_protti, eval=FALSE----------------------------------------------
#  # Load data
#  your_data <- read_protti("your_differential_abundance_data.csv")

## ----prepare_data, eval = test_protti, warning = FALSE------------------------
#  # Input UniProt IDs
#  uniprot_ids <- unique(ptsi_pgk$pg_protein_accessions)
#  
#  # Fetch UniProt information
#  uniprot_information <- fetch_uniprot(uniprot_ids = uniprot_ids,
#                                       columns = c("sequence", "database(PDB)"))
#  
#  # Add UniProt information and find peptide positions
#  ptsi_pgk_annotated <- ptsi_pgk %>%
#    left_join(uniprot_information, by = c("pg_protein_accessions" = "id")) %>%
#    find_peptide(protein_sequence = sequence, peptide_sequence = pep_stripped_sequence)

## ----extract_pdb_info, eval = test_protti, warning = FALSE--------------------
#  # Extract PDB IDs from UniProt information
#  ptsi_pgk_pdb_ids <- ptsi_pgk_annotated %>%
#    distinct(pg_protein_accessions, database_pdb) %>%
#    mutate(pdb_id = str_split(database_pdb, pattern = ";")) %>%
#    unnest(pdb_id) %>%
#    filter(pdb_id != "")
#  
#  # Fetch pdb information
#  ptsi_pgk_pdb_information <- fetch_pdb(pdb_ids = unique(ptsi_pgk_pdb_ids$pdb_id))

## ----filter_structures, eval = test_protti, warning = FALSE-------------------
#  filtered_structures <- ptsi_pgk_pdb_information %>%
#    filter(experimental_method == "X-ray",
#           resolution_combined <= 3) %>%
#    group_by(reference_database_accession) %>%
#    filter(length == max(length)) %>%
#    ungroup()

## ----fetch_pdb_structure, eval = test_protti, warning = FALSE-----------------
#  # Input PDB IDs
#  pdb_ids <- unique(filtered_structures$pdb_ids) # "1ZMR", "2HWG"
#  
#  # Fetch atom level structural information
#  ptsi_pgk_structure_information <- fetch_pdb_structure(pdb_ids = pdb_ids,
#                                                        return_data_frame = TRUE)

## ----fetch_alphafold_prediction, eval = test_protti, warning = FALSE----------
#  # Fetch atom level structural prediction information from AlphaFold
#  ptsi_pgk_prediction_information <- fetch_alphafold_prediction(uniprot_ids = uniprot_ids,
#                                                                return_data_frame = TRUE)
#  
#  # Example for fetching all predictions for Methanocaldococcus jannaschii
#  # mj_predictions <- fetch_alphafold_prediction(organism_name = "Methanocaldococcus jannaschii")

## ----find_peptide, eval = test_protti, warning = FALSE------------------------
#  ptsi_pgk_peptide_structure_positions <- find_peptide_in_structure(
#     peptide_data = ptsi_pgk_annotated,
#     peptide = pep_stripped_sequence,
#     start = start,
#     end = end,
#     uniprot_id = pg_protein_accessions,
#     pdb_data = filtered_structures,
#     retain_columns = c(eg_precursor_id, diff, adj_pval))

## ----create_structure_contact_map, eval = test_protti, warning = FALSE, fig.width = 10, fig.height = 7, fig.align = "center"----
#  # Filter data for significant peptides.
#  significant_peptides <- ptsi_pgk_peptide_structure_positions %>%
#    filter(abs(diff) > 2, adj_pval <= 0.01)
#  
#  # Create a structure contact maps
#  contact_map <- create_structure_contact_map(
#    data = significant_peptides,
#    id = pdb_ids,
#    chain = auth_asym_id,
#    auth_seq_id = auth_seq_id,
#    distance_cutoff = 10,
#    pdb_model_number_selection = c(0, 1),
#    return_min_residue_distance = TRUE
#  )
#  
#  # This is a helper function for the plot.
#  # It allows the display of integers on the axis.
#  integer_breaks <- function(n = 5, ...) {
#    fxn <- function(x) {
#      breaks <- floor(pretty(x, n, ...))
#      names(breaks) <- attr(breaks, "labels")
#      breaks
#    }
#    return(fxn)
#  }
#  
#  # Plot structure contact maps
#  # 1ZMR
#  contact_map[["1ZMR"]] %>% # Extract data frame from list
#    mutate(chain_combinations = paste0("chain_", label_asym_id_var1, "_vs_chain_", label_asym_id_var2)) %>%
#    ggplot(aes(x = label_seq_id_var1, y = label_seq_id_var2, fill = min_distance_residue)) +
#    geom_tile() +
#    scale_y_continuous(breaks = integer_breaks()) +
#    scale_x_continuous(breaks = integer_breaks()) +
#    facet_wrap(~chain_combinations, scale = "free") +
#    labs(title = "Structure contact map 1ZMR") +
#    theme_bw()
#  
#  # 2HWG
#  contact_map[["2HWG"]] %>% # Extract data frame from list
#    mutate(chain_combinations = paste0("chain_", label_asym_id_var1, "_vs_chain_", label_asym_id_var2)) %>%
#    ggplot(aes(x = label_seq_id_var1, y = label_seq_id_var2, fill = min_distance_residue)) +
#    geom_tile() +
#    scale_y_continuous(breaks = integer_breaks()) +
#    scale_x_continuous(breaks = integer_breaks()) +
#    labs(title = "Structure contact map 2HWG") +
#    facet_wrap(~chain_combinations, scale = "free") +
#    theme_bw()

## ----peptide_mapping, eval = test_protti, warning = FALSE---------------------
#  ptsi_pgk_peptide_structure_positions %>%
#    mutate(map_value = ifelse(eg_precursor_id %in% significant_peptides$eg_precursor_id,
#                              100,
#                              0)) %>%
#    map_peptides_on_structure(
#     uniprot_id = pg_protein_accessions,
#     pdb_id = pdb_ids,
#     chain = auth_asym_id,
#     auth_seq_id = auth_seq_id,
#     map_value = map_value,
#     file_format = ".pdb",
#     export_location = tempdir() # change to a location of your choice
#   )

## ----3d_structure_mapping, eval = test_protti, echo=TRUE, warning=FALSE-------
#  # Install the r3dmol package if it is not installed
#  # install.packages("r3dmol")
#  
#  # Load the r3dmol package
#  library(r3dmol)
#  
#  # Create structure
#  r3dmol() %>%
#    m_add_model(data = paste0(tempdir(), "/1ZMR_P0A799.pdb"), format = "pdb") %>%
#    m_set_style(style = m_style_cartoon(
#      colorfunc = "
#          function(atom) {
#            if (atom.b == 50) {return '#90EE90'};
#            if (atom.b == 100) {return '#FF7276'};
#            if (atom.b == 0) {return 'white'};
#            return 'white';
#          }"
#    )) %>%
#    m_zoom_to()

## ----peptide_map_2hwg, eval = test_protti, echo = FALSE, fig.align = "center", fig.cap = "**ptsI (2HWG) with mapped significantly changing peptides.**"----
#  knitr::include_graphics("figures/peptide_map_2hwg.png")

## ----interaction_2hwg, eval = test_protti, echo = FALSE, fig.align = "center", out.width = "60%", fig.cap = "**ptsI (2HWG) binding interface of significantly chaging peptide.**"----
#  knitr::include_graphics("figures/interaction_2hwg.png")

## ----peptide_map_1zmr, eval = test_protti, echo = FALSE, fig.align = "center", fig.cap = "**pgk (1ZMR) with mapped significantly changing peptides.**"----
#  knitr::include_graphics("figures/peptide_map_1zmr.png")

## ----calculate_and_map_scores, eval = test_protti-----------------------------
#  # Calculate the amino acid score
#  amino_acid_score <- calculate_aa_scores(
#    data = ptsi_pgk_peptide_structure_positions,
#    protein = pg_protein_accessions,
#    diff = diff,
#    adj_pval = adj_pval,
#    start_position = start,
#    end_position = end,
#    retain_columns = c(pdb_ids, auth_asym_id)
#  )
#  
#  # Find amino acid positions in the structure
#  ptsi_pgk_amino_acid_structure_positions <- find_peptide_in_structure(
#     peptide_data = amino_acid_score,
#     peptide = residue,
#     start = residue,
#     end = residue,
#     uniprot_id = pg_protein_accessions,
#     pdb_data = filtered_structures,
#     retain_columns = c(amino_acid_score))
#  
#  # Map the score on structure
#  map_peptides_on_structure(
#    peptide_data = ptsi_pgk_amino_acid_structure_positions,
#     uniprot_id = pg_protein_accessions,
#     pdb_id = pdb_ids,
#     chain = auth_asym_id,
#     auth_seq_id = auth_seq_id,
#     map_value = amino_acid_score,
#     file_format = ".pdb",
#     export_location = tempdir()
#   )

## ----score_3d_structure_mapping, eval = test_protti, echo=TRUE, warning=FALSE----
#  # create a color gradient with 101 colors
#  color_gradient <- paste0('"',
#                           paste(colorRampPalette(c("white", "#90EE90", "#FF7276"))(101),
#                                 collapse = '", "'),
#                           '"')
#  
#  # create structure
#  r3dmol() %>%
#    m_add_model(data = paste0(tempdir(), "/2HWG_P08839.pdb"), format = "pdb") %>%
#    m_set_style(style = m_style_cartoon(
#      colorfunc = paste0("
#          function(atom) {
#            const color = [", color_gradient,"]
#            return color[Math.round(atom.b)]
#          }")
#    )) %>%
#    m_zoom_to()

## ----peptide_map_2hwg_score, eval = test_protti, echo = FALSE, fig.align = "center", fig.cap = "**ptsI (2HWG) with mapped amino acid scores.**"----
#  knitr::include_graphics("figures/peptide_map_2hwg_score.png")

## ----peptide_map_1zmr_score, eval = test_protti, echo = FALSE, fig.align = "center", fig.cap = "**pgk (1ZMR) with mapped amino acid scores.**"----
#  knitr::include_graphics("figures/peptide_map_1zmr_score.png")

