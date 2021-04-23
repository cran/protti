## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)

## ----Spectronaut, eval=FALSE--------------------------------------------------
#  # To read in your own data you can use read_protti()
#  spectronaut_data  <- read_protti(filename = "mydata/spectronaut.csv")

## ----MaxQuant_peptide, eval=FALSE---------------------------------------------
#  # To read in your own data you can use read_protti()
#  evidence  <- read_protti(filename = "yourpath/evidence.txt")
#  
#  evidence_proteotypic <- evidence %>%
#    mutate(is_proteotypic = str_detect(
#      string = proteins,
#      pattern = ";",
#      negate = TRUE
#    )) %>%# adds new column with logicals that are TRUE if the peptide can be assigned to only one protein and FALSE if it can be assigned to multiple
#    mutate(is_contaminant = ifelse(potential_contaminant == "+", TRUE, FALSE)) # adds new column with logicals indicating if peptide is coming from a potential contaminant
#  
#  # Make an annotation data frame and merge it with your data frame to obtain conditions
#  # We are annotating sample 1-3 as controls and samples 4-6 as treated conditions
#  
#  file_name <- c( # make sure that the names are the same name as in your report
#    "sample1",
#    "sample2",
#    "sample3",
#    "sample4",
#    "sample5",
#    "sample6"
#  )
#  
#  condition <- c(
#    "control",
#    "control",
#    "control",
#    "treated",
#    "treated",
#    "treated"
#  )
#  
#  annotation <- data.frame(file_name, condition)
#  
#  # Combine your long data frame with the annotation
#  evidence_annotated <- evidence_proteotypic %>%
#    left_join(y = annotation, by = "file_name")

## ----MaxQuant_protein, eval=FALSE---------------------------------------------
#  # To read in your own data you can use read_protti()
#  protein_groups  <- read_protti(filename = "yourpath/proteinGroups.txt") %>%
#    mutate(is_potential_contaminant = ifelse(potential_contaminant == "+", TRUE, FALSE)) # adds new column with logicals indicating if protein is a potential contaminant, you can filter these out later on. You should also consider filtering out proteins that were "only identified by site" and reverse hits, as well as proteins with only one identified peptide
#  
#  # Change wide format to long format and create new columns called `r_file_name`and `intensity`
#  protein_groups_long <- protein_groups %>%
#    pivot_longer(
#      cols = starts_with("intensity_"),
#      names_to = "file_name",
#      values_to = "intensity"
#    )
#  
#  # Make an annotation data frame and merge it with your data frame to obtain conditions
#  # We are annotating sample 1-3 as controls and samples 4-6 as treated conditions
#  
#  file_name <- c( # make sure that the names are the same name as in your report
#    "intensity_sample1",
#    "intensity_sample2",
#    "intensity_sample3",
#    "intensity_sample4",
#    "intensity_sample5",
#    "intensity_sample6"
#  )
#  
#  condition <- c(
#    "control",
#    "control",
#    "control",
#    "treated",
#    "treated",
#    "treated"
#  )
#  
#  annotation <- data.frame(file_name, condition)
#  
#  # Combine your long data frame with the annotation
#  protein_groups_annotated <- protein_groups_long %>%
#    left_join(y = annotation, by = "file_name")

## ----Skyline, eval=FALSE------------------------------------------------------
#  # Load data
#  skyline_data <- read_protti(filename = "yourpath/skyline.csv")
#  
#  skyline_data_int <- skyline_data %>%
#    mutate(precursor = paste0(peptide_sequence, "_", charge)) %>% # create a column with precursor information
#    group_by(replicate_name, precursor) %>%
#    mutate(sum_intensity = sum(area)) %>% # making a new column containing the summed up intensities of all transitions of one precursor
#    select(-c(product_mz, area)) %>% # removing the columns we don't need
#    distinct() # removing duplicated rows from the data frame
#  
#  # Add annotation
#  replicate_name <- c( # make sure that the names are the same name as in your report
#    "sample_1",
#    "sample_2",
#    "sample_3",
#    "sample_1",
#    "sample_2",
#    "sample_3"
#  )
#  
#  condition <- c(
#    "control",
#    "control",
#    "control",
#    "treated",
#    "treated",
#    "treated"
#  )
#  
#  annotation <- data.frame(replicate_name, condition)
#  
#  # Combine your long data frame with the annotation
#  skyline_annotated <- skyline_data_int %>%
#    left_join(y = annotation, by = "replicate_name")
#  

## ----Proteome_discoverer_pep, eval=FALSE--------------------------------------
#  # Load data
#  pd_pep_data <- read_protti("yourpath/PDpeptides.csv")
#  
#  # Select relevant columns
#  pd_pep_selected <- pd_pep_data %>%
#    select(
#      sequence,
#      modifications,
#      number_proteins,
#      contaminant,
#      master_protein_accessions,
#      starts_with("abundances_grouped"), # select all columns that start with "abundances_grouped"
#      quan_info
#    )
#  
#  # Filter data frame
#  pd_pep_filtered <- pd_pep_selected %>%
#    filter(contaminant == FALSE) %>% # remove annotated contaminants
#    filter(number_proteins == 1) %>% # select proteotypic peptides
#    filter(quan_info != "No Quan Values") # remove peptides that have no quantification values
#  
#  # Convert into long format
#  pd_pep_long <- pd_pep_filtered %>%
#    pivot_longer(
#      cols = starts_with("abundances"),
#      names_to = "file_name",
#      values_to = "intensity"
#    ) %>%
#    mutate(precursor = paste(sequence, modifications)) # combine peptide sequence and modifications to make a precursor column
#  
#  # Make annotation data frame
#  file_name <- c( # make sure that the names are the same name as in your report
#    "abundances_grouped_f1",
#    "abundances_grouped_f2",
#    "abundances_grouped_f3",
#    "abundances_grouped_f4",
#    "abundances_grouped_f5",
#    "abundances_grouped_f6"
#  )
#  
#  condition <- c(
#    "control",
#    "control",
#    "control",
#    "treated",
#    "treated",
#    "treated"
#  )
#  
#  annotation <- data.frame(file_name, condition)
#  
#  # Combine your long data frame with the annotation
#  pd_pep_long_annotated <- pd_pep_long %>%
#    left_join(y = annotation, by = "file_name")

## ----Proteome_discoverer_prot, eval=FALSE-------------------------------------
#  # Load data
#  pd_prot_data <- read_protti("yourpath/PDproteins.csv")
#  
#  # Select relevant columns
#  pd_prot_selected <- pd_prot_data %>%
#    select(
#      accession,
#      description,
#      contaminant,
#      number_peptides,
#      starts_with("abundances_grouped"), # select all columns that start with "abundances_grouped"
#    )
#  
#  # Filter data frame
#  pd_prot_data_filtered <- pd_prot_selected %>%
#    filter(contaminant == FALSE) %>% # remove annotated contaminants
#    filter(number_peptides > 1) # select proteins with more than one identified peptide
#  
#  # Convert into long format
#  pd_prot_long <- pd_prot_data_filtered %>%
#    pivot_longer(
#      cols = starts_with("abundances"),
#      names_to = "file_name",
#      values_to = "intensity"
#    )
#  
#  # Make annotation data frame
#  file_name <- c( # make sure that the names are the same name as in your report
#    "abundances_grouped_f1",
#    "abundances_grouped_f2",
#    "abundances_grouped_f3",
#    "abundances_grouped_f4",
#    "abundances_grouped_f5",
#    "abundances_grouped_f6"
#  )
#  
#  condition <- c(
#    "control",
#    "control",
#    "control",
#    "treated",
#    "treated",
#    "treated"
#  )
#  
#  annotation <- data.frame(file_name, condition)
#  
#  # Combine your long data frame with the annotation
#  pd_prot_long_annotated <- pd_prot_long %>%
#    left_join(y = annotation, by = "file_name")

