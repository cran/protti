% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomise_queue.R
\name{randomise_queue}
\alias{randomise_queue}
\title{Randomise samples in MS queue}
\usage{
randomise_queue(data = NULL, rows = NULL, export = FALSE)
}
\arguments{
\item{data}{optional, a data frame that contains a queue. If not provided a queue file can be
chosen interactively.}

\item{rows}{optional, a numeric vector that specifies a range of rows in for which samples
should be randomized.}

\item{export}{a logical value that determines if a \code{"randomised_queue.csv"} file will be
saved in the working directory. If FALSE a data frame will be returned.}
}
\value{
If \code{export = TRUE} a \code{"randomised_queue.csv"} file will be saved in the
working directory. If \code{export = FALSE} a data frame that contains the randomised queue
is returned.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
This function randomises the order of samples in an MS queue. QC and Blank samples are left in
place. It is also possible to randomise only parts of the queue. Before running this make sure
to set a specific seed with the \code{set.seed()} function. This ensures that the randomisation
of the result is consistent if the function is run again.
}
\examples{
queue <- create_queue(
  date = c("200722"),
  instrument = c("EX1"),
  user = c("jquast"),
  measurement_type = c("DIA"),
  experiment_name = c("JPQ031"),
  digestion = c("LiP", "tryptic control"),
  treatment_type_1 = c("EDTA", "H2O"),
  treatment_type_2 = c("Zeba", "unfiltered"),
  treatment_dose_1 = c(10, 30, 60),
  treatment_unit_1 = c("min"),
  n_replicates = 4,
  number_runs = FALSE,
  organism = c("E. coli"),
  exclude_combinations = list(list(
    treatment_type_1 = c("H2O"),
    treatment_type_2 = c("Zeba", "unfiltered"),
    treatment_dose_1 = c(10, 30)
  )),
  inj_vol = c(2),
  data_path = "D:\\\\2007_Data",
  method_path = "C:\\\\Xcalibur\\\\methods\\\\DIA_120min",
  position_row = c("A", "B", "C", "D", "E", "F"),
  position_column = 8,
  blank_every_n = 4,
  blank_position = "1-V1",
  blank_method_path = "C:\\\\Xcalibur\\\\methods\\\\blank"
)

head(queue, n = 20)

randomised_queue <- randomise_queue(
  data = queue,
  export = FALSE
)

head(randomised_queue, n = 20)
}
