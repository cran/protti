#' Fetch PDB structure atom data from RCSB
#'
#' Fetches atom data for a PDB structure from RCSB. If you want to retrieve metadata about PDB
#' structures, use the function \code{fetch_pdb()}. The information retrieved is based on the
#' .cif file of the structure, which may vary from the .pdb file.
#'
#' @param pdb_ids a character vector of PDB identifiers.
#' @param return_data_frame a logical value that indicates if a data frame instead of a list is
#' returned. It is recommended to only use this if not many pdb structures are retrieved. Default
#' is FALSE.
#' @param show_progress a logical value that indicates if a progress bar will be shown.
#' Default is TRUE.
#'
#' @return A list that contains atom data for each PDB structures provided. If return_data_frame is
#' TRUE, a data frame with this information is returned instead. The data frame contains the
#' following columns:
#'
#' * label_id: Uniquely identifies every atom in the structure following the standardised
#' convention for mmCIF files. Example value: "5", "C12", "Ca3g28", "Fe3+17", "H*251", "boron2a",
#' "C a phe 83 a 0", "Zn Zn 301 A 0"
#' * type_symbol: The code used to identify the atom species representing this atom type.
#' Normally this code is the element symbol. The code may be composed of any character except an
#' underscore with the additional proviso that digits designate an oxidation state and must be
#' followed by a + or - character. Example values: "C", "Cu2+", "H(SDS)", "dummy", "FeNi".
#' * label_atom_id: Uniquely identifies every atom for the given residue following the
#' standardised convention for mmCIF files. Example values: "CA", "HB1", "CB", "N"
#' * label_comp_id: A chemical identifier for the residue. For protein polymer entities,
#' this is the three- letter code for the amino acid. For nucleic acid polymer entities, this is
#' the one-letter code for the base. Example values: "ala", "val", "A", "C".
#' * label_asym_id: Chain identifier following the standardised convention for mmCIF files.
#' Example values: "1", "A", "2B3".
#' * entity_id: Records details about the molecular entities that are present in the
#' crystallographic structure. Usually all different types of molecular entities such as polymer
#' entities, non-polymer entities or water molecules are numbered once for each structure. Each
#' type of non-polymer entity has its own number. Thus, the highest number in this column
#' represents the number of different molecule types in the structure.
#' * label_seq_id: Uniquely and sequentially identifies residues for each \code{label_asym_id}.
#' This is always a number and the sequence of numbers always progresses in increasing numerical order.
#' * x: The x coordinate of the atom.
#' * y: The y coordinate of the atom.
#' * z: The z coordinate of the atom.
#' * site_occupancy: The fraction of the atom type present at this site.
#' * b_iso_or_equivalent: Contains the B-factor or isotopic atomic displacement factor for
#' each atom.
#' * formal_charge: The net integer charge assigned to this atom. This is the formal charge
#' assignment normally found in chemical diagrams. It is currently only assigned in a small subset
#' of structures.
#' * auth_seq_id: An alternative residue identifier (\code{label_seq_id}) provided by the
#' author of the structure in order to match the identification used in the publication that
#' describes the structure. This does not need to be numeric and is therefore of type character.
#' * auth_comp_id: An alternative chemical identifier (\code{label_comp_id}) provided by the
#' author of the structure in order to match the identification used in the publication that
#' describes the structure.
#' * auth_asym_id: An alternative chain identifier (\code{label_asym_id}) provided by the
#' author of the structure in order to match the identification used in the publication that
#' describes the structure.
#' * pdb_model_number: The PDB model number.
#' * pdb_id: The protein database identifier for the structure.
#'
#'
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom stringr str_replace_all str_detect
#' @importFrom curl has_internet
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @export
#'
#' @examples
#' \donttest{
#' pdb_structure <- fetch_pdb_structure(
#'   pdb_ids = c("6HG1", "1E9I", "6D3Q", "4JHW"),
#'   return_data_frame = TRUE
#' )
#'
#' head(pdb_structure, n = 10)
#' }
fetch_pdb_structure <- function(pdb_ids, return_data_frame = FALSE, show_progress = TRUE) {
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  # remove NA values
  pdb_ids <- pdb_ids[!is.na(pdb_ids)]

  batches <- purrr::map(
    .x = pdb_ids,
    .f = ~ paste0("https://files.rcsb.org/download/", .x, ".cif")
  )

  names(batches) <- pdb_ids

  if (show_progress == TRUE) {
    pb <- progress::progress_bar$new(
      total = length(batches),
      format = "  Fetching structures [:bar] :current/:total (:percent) :eta"
    )
  }

  query_result <- purrr::map2(
    .x = batches,
    .y = names(batches),
    .f = ~ {
      # query information from database
      query <- try_query(.x,
        type = "text/tab-separated-values",
        col_names = FALSE,
        quote = "",
        show_col_types = FALSE,
        progress = FALSE
      )

      if (show_progress == TRUE) {
        pb$tick()
      }

      # only proceed with data if it was correctly retrieved
      if ("tbl" %in% class(query)) {
        query %>%
          dplyr::filter(stringr::str_detect(X1,
            pattern = "^ATOM\\s+\\d|^HETATM\\s+\\d"
          )) %>%
          dplyr::mutate(X2 = stringr::str_replace_all(X1,
            pattern = "\\s+",
            replacement = " "
          )) %>%
          tidyr::separate(X2,
            sep = " ",
            into = c(
              "x1",
              "label_id",
              "type_symbol",
              "label_atom_id",
              "x2",
              "label_comp_id",
              "label_asym_id",
              "entity_id",
              "label_seq_id",
              "x3",
              "x",
              "y",
              "z",
              "site_occupancy",
              "b_iso_or_equivalent",
              "formal_charge",
              "auth_seq_id",
              "auth_comp_id",
              "auth_asym_id",
              "x4",
              "pdb_model_number"
            )
          ) %>%
          dplyr::select(-c("X1", "x1", "x2", "x3", "x4")) %>%
          dplyr::group_by(.data$label_asym_id, .data$label_atom_id, .data$label_comp_id) %>%
          dplyr::mutate(label_seq_id = ifelse(.data$label_seq_id == ".",
            1:n(),
            as.numeric(.data$label_seq_id)
          )) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            label_id = as.numeric(.data$label_id),
            entity_id = as.numeric(.data$entity_id),
            label_seq_id = as.numeric(.data$label_seq_id),
            x = as.numeric(.data$x),
            y = as.numeric(.data$y),
            z = as.numeric(.data$z),
            site_occupancy = as.numeric(.data$site_occupancy),
            b_iso_or_equivalent = as.numeric(.data$b_iso_or_equivalent),
            auth_seq_id = .data$auth_seq_id,
            pdb_model_number = as.numeric(.data$pdb_model_number),
            pdb_id = .y
          )
      } else {
        query
      }
    }
  )

  # catch any IDs that have not been fetched correctly
  error_list <- query_result %>%
    purrr::keep(.p = ~ is.character(.x))

  error_table <- tibble::tibble(
    id = names(error_list),
    error = unlist(error_list)
  ) %>%
    dplyr::distinct()

  if (nrow(error_table) != 0) {
    message("The following IDs have not be retrieved correctly.")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # only keep data in output

  query_result <- query_result %>%
    purrr::keep(.p = ~ !is.character(.x))

  if (return_data_frame == FALSE) {
    return(query_result)
  } else {
    query_result_df <- purrr::map_dfr(
      .x = query_result,
      .f = ~.x
    )
    return(query_result_df)
  }
}
