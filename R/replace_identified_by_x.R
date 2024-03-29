#' Replace identified positions in protein sequence by "x"
#'
#' Helper function for the calculation of sequence coverage, replaces identified positions with an
#' "x" within the protein sequence.
#'
#' @param sequence a character value that contains the protein sequence.
#' @param positions_start a numeric vector of start positions of the identified peptides.
#' @param positions_end a numeric vector of end positions of the identified peptides.
#'
#' @return A character vector that contains the modified protein sequence with each identified
#' position replaced by "x".
#' @importFrom purrr map2
#' @importFrom stringr str_sub
replace_identified_by_x <-
  function(sequence, positions_start, positions_end) {
    sequence <- unique(sequence)
    if (sequence == "" | is.na(sequence)) {
      return(NA)
    }
    remove_na <- !is.na(positions_start) & !is.na(positions_end)
    positions_start <- positions_start[remove_na]
    positions_end <- positions_end[remove_na]
    result <- purrr::map2(
      .x = positions_start, .y = positions_end,
      function(x, y) {
        times <- y - x + 1
        stringr::str_sub(sequence, start = x, end = y) <- paste(rep("x", times = times), collapse = "")
        # this does not modify the global environment but only the
        # environment of the parent function (replace_identified_by_x).
        sequence <<- sequence
      }
    )
    result[[length(result)]]
  }
