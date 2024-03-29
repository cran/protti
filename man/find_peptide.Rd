% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_peptide.R
\name{find_peptide}
\alias{find_peptide}
\title{Find peptide location}
\usage{
find_peptide(data, protein_sequence, peptide_sequence)
}
\arguments{
\item{data}{a data frame that contains at least the protein and peptide sequence.}

\item{protein_sequence}{a character column in the \code{data} data frame that contains the
protein sequence.}

\item{peptide_sequence}{a character column in the \code{data} data frame that contains the
peptide sequence.}
}
\value{
A data frame that contains the input data and four additional columns with peptide
start and end position, the last amino acid and the amino acid before the peptide.
}
\description{
The position of the given peptide sequence is searched within the given protein sequence. In
addition the last amino acid of the peptide and the amino acid right before are reported.
}
\examples{
# Create example data
data <- data.frame(
  protein_sequence = c("abcdefg"),
  peptide_sequence = c("cde")
)

# Find peptide
find_peptide(
  data = data,
  protein_sequence = protein_sequence,
  peptide_sequence = peptide_sequence
)
}
