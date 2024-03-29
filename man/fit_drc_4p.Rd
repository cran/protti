% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_drc_4p.R
\name{fit_drc_4p}
\alias{fit_drc_4p}
\title{Fitting four-parameter dose response curves}
\usage{
fit_drc_4p(
  data,
  sample,
  grouping,
  response,
  dose,
  filter = "post",
  replicate_completeness = 0.7,
  condition_completeness = 0.5,
  correlation_cutoff = 0.8,
  log_logarithmic = TRUE,
  include_models = FALSE,
  retain_columns = NULL
)
}
\arguments{
\item{data}{a data frame that contains at least the input variables.}

\item{sample}{a character column in the \code{data} data frame that contains the sample names.}

\item{grouping}{a character column in the \code{data} data frame that contains the precursor,
peptide or protein identifiers.}

\item{response}{a numeric column in the \code{data} data frame that contains the response
values, e.g. log2 transformed intensities.}

\item{dose}{a numeric column in the \code{data} data frame that contains the dose values, e.g.
the treatment concentrations.}

\item{filter}{a character value that determines if models should be filtered and if they should
be filtered before or after the curve fits. Filtering of models can be skipped with
\code{filter = "none"}. Data can be filtered prior to model fitting with \code{filter = "pre"}.
In that case models will only be fitted for data that passed the filtering step. This will
allow for faster model fitting since only fewer models will be fit. If you plan on performing
an enrichment analysis you have to choose \code{filter = "post"}. All models will be fit (even
the ones that do not pass the filtering criteria). For enrichment analysis you should use both
good (i.e. models that pass the filtering) and bad (i.e. models that do not pass the filtering)
models. Therefore, for post-filtering the full list is returned and it will only contain
annotations that indicate (\code{passed_filter}) if the filtering was passed or not. Default is
"post". For ANOVA an adjusted p-value of 0.05 is used as a cutoff.}

\item{replicate_completeness}{a numeric value which similar to \code{completenss_MAR} of the
\code{assign_missingness} function sets a threshold for the completeness of data. In contrast
to \code{assign_missingness} it only determines the completeness for one condition and not the
comparison of two conditions. The threshold is used to calculate a minimal degree of data
completeness. The value provided to this argument has to be between 0 and 1, default is 0.7.
It is multiplied with the number of replicates and then adjusted downward. The resulting number
is the minimal number of observations that a condition needs to have to be considered "complete
enough" for the \code{condition_completeness} argument.}

\item{condition_completeness}{a numeric value which determines how many conditions need to at
least fulfill the "complete enough" criteria set with \code{replicate_completeness}. The
value provided to this argument has to be between 0 and 1, default is 0.5. It is multiplied with
the number of conditions and then adjusted downward. The resulting number is the minimal number
of conditions that need to fulfill the \code{replicate_completeness} argument for a peptide to
pass the filtering.}

\item{correlation_cutoff}{a numeric vector that specifies the correlation cutoff used for data
filtering.}

\item{log_logarithmic}{a logical value that indicates if a logarithmic or log-logarithmic model
is fitted. If response values form a symmetric curve for non-log transformed dose values, a
logarithmic model instead of a log-logarithmic model should be used. Usually biological dose
response data has a log-logarithmic distribution, which is the reason this is the default.
Log-logarithmic models are symmetric if dose values are log transformed.}

\item{include_models}{a logical value that indicates if model fit objects should be exported.
These are usually very large and not necessary for further analysis.}

\item{retain_columns}{a vector that specifies columns that should be retained from the input
data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
columns can be retained by providing their names (not in quotations marks, just like other
column names, but in a vector).}
}
\value{
If \code{include_models = FALSE} a data frame is returned that contains correlations
of predicted to measured values as a measure of the goodness of the curve fit, an associated
p-value and the four parameters of the model for each group. Furthermore, input data for plots
is returned in the columns \code{plot_curve} (curve and confidence interval) and \code{plot_points}
(measured points). If \code{include_models = TURE}, a list is returned that contains:
\itemize{
\item \code{fit_objects}: The fit objects of type \code{drc} for each group.
\item \code{correlations}: The correlation data frame described above
}
}
\description{
Function for fitting four-parameter dose response curves for each group (precursor, peptide or
protein). In addition it can filter data based on completeness, the completeness distribution
and statistical testing using ANOVA.
}
\details{
If data filtering options are selected, data is filtered based on multiple criteria.
In general, curves are only fitted if there are at least 5 conditions with data points present
to ensure that there is potential for a good curve fit. Therefore, this is also the case if no
filtering option is selected. Furthermore, a completeness cutoff is defined for filtering. By
default each entity (e.g. precursor) is filtered to contain at least 70\% of total replicates
(adjusted downward) for at least 50\% of all conditions (adjusted downward). This can be adjusted
with the according arguments. In addition to the completeness cutoff, also a significance cutoff
is applied. ANOVA is used to compute the statistical significance of the change for each entity.
The resulting p-value is adjusted using the Benjamini-Hochberg method and a cutoff of q <= 0.05
is applied. Curve fits that have a minimal value that is higher than the maximal value are
excluded as they were likely wrongly fitted. Curves with a correlation below 0.8 are not passing
the filtering. If a fit does not fulfill the significance or completeness cutoff, it has a chance
to still be considered if half of its values (+/-1 value) pass the replicate completeness
criteria and half do not pass it. In order to fall into this category, the values that fulfill t
he completeness cutoff and the ones that do not fulfill it need to be consecutive, meaning
located next to each other based on their concentration values. Furthermore, the values that
do not pass the completeness cutoff need to be lower in intensity. Lastly, the difference
between the two groups is tested for statistical significance using a Welch's t-test and a
cutoff of p <= 0.1 (we want to mainly discard curves that falsly fit the other criteria but that
have clearly non-significant differences in mean). This allows curves to be considered that have
missing values in half of their observations due to a decrease in intensity. It can be thought
of as conditions that are missing not at random (MNAR). It is often the case that those entities
do not have a significant p-value since half of their conditions are not considered due to data
missingness.

The final filtered list is ranked based on a score calculated on entities that pass the filter.
The score is the negative log10 of the adjusted ANOVA p-value scaled between 0 and 1 and the
correlation scaled between 0 and 1 summed up and divided by 2. Thus, the highest score an
entity can have is 1 with both the highest correlation and adjusted p-value. The rank is
corresponding to this score. Please note, that entities with MNAR conditions might have a
lower score due to the missing or non-significant ANOVA p-value. You should have a look at
curves that are TRUE for \code{dose_MNAR} in more detail.
}
\examples{
\donttest{
# Load libraries
library(dplyr)

set.seed(123) # Makes example reproducible

# Create example data
data <- create_synthetic_data(
  n_proteins = 2,
  frac_change = 1,
  n_replicates = 3,
  n_conditions = 8,
  method = "dose_response",
  concentrations = c(0, 1, 10, 50, 100, 500, 1000, 5000),
  additional_metadata = FALSE
)

# Perform dose response curve fit
drc_fit <- fit_drc_4p(
  data = data,
  sample = sample,
  grouping = peptide,
  response = peptide_intensity_missing,
  dose = concentration,
  retain_columns = c(protein, change_peptide)
)

glimpse(drc_fit)

head(drc_fit, n = 10)
}
}
