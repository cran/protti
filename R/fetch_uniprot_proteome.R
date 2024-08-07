#' Fetch proteome data from UniProt
#'
#' Fetches proteome data from UniProt for the provided organism ID.
#'
#' @param organism_id a numeric value that specifies the NCBI taxonomy identifier (TaxId) for an
#' organism.
#' @param columns a character vector of metadata columns that should be imported from UniProt (all
#' possible columns can be found \href{https://www.uniprot.org/help/return_fields}{here}. For
#' cross-referenced database provide the database name with the prefix "xref_", e.g. \code{"xref_pdb"}).
#' Note: Not more than one or two columns should be selected otherwise the function will not be
#' able to efficiently retrieve the information. If more information is needed, \code{fetch_uniprot()}
#' can be used with the IDs retrieved by this function.
#' @param reviewed a logical value that determines if only reviewed protein entries will be retrieved.
#' @param timeout a numeric value specifying the time in seconds until the download times out.
#' The default is 60 seconds.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. The default is 2.
#'
#' @return A data frame that contains all protein metadata specified in \code{columns} for the
#' organism of choice.
#' @importFrom janitor make_clean_names
#' @export
#'
#' @examples
#' \donttest{
#' head(fetch_uniprot_proteome(9606))
#' }
fetch_uniprot_proteome <-
  function(organism_id,
           columns = c("accession"),
           reviewed = TRUE,
           timeout = 120,
           max_tries = 5) {
    if (!curl::has_internet()) {
      message("No internet connection.")
      return(invisible(NULL))
    }

    if (length(organism_id) == 0) {
      stop("No valid organism ID found.")
    }
    if (length(columns) > 4) {
      warning(strwrap("We suggest to use the fetch_uniprot function to fetch more than four columns.",
        prefix = "\n", initial = ""
      ))
    }
    url <- "http://rest.uniprot.org/uniprotkb/stream?query="
    column_names <- janitor::make_clean_names(columns)
    collapsed_columns <- paste(columns, collapse = ",")
    reviewed <- paste0("reviewed:", ifelse(reviewed == TRUE, "true", "false"))
    organism_id <- paste0("organism_id:", organism_id)
    query_url <-
      utils::URLencode(paste0(
        url,
        reviewed,
        "+AND+",
        organism_id,
        "&format=tsv&fields=",
        collapsed_columns
      ))
    result <- try_query(query_url, timeout = timeout, max_tries = max_tries, silent = FALSE, progress = FALSE, show_col_types = FALSE)
    # result can either be a data.frame or it is a character string with the error message
    if (!methods::is(result, "data.frame")) {
      if (stringr::str_detect(result, pattern = "Timeout")) {
        message('The data retrieval timed out. Consider increasing the "timeout" or "max_tries" argument. \n')
      }
      return(invisible(result))
    }
    colnames(result) <- column_names
    result
  }
