% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregateByAttributes.R
\name{aggregateByAttributes}
\alias{aggregateByAttributes}
\title{Aggregate PSMs to Protein Level}
\usage{
aggregateByAttributes(
  PSM,
  colOfReporterIonInt,
  ranger = NULL,
  predError = NULL,
  ratioCalc = FALSE
)
}
\arguments{
\item{PSM}{A data frame containing all PSMs to be aggregated.}

\item{colOfReporterIonInt}{A vector of column names representing reporter ion intensities across different channels.}

\item{ranger}{The random forest model to be applied for aggregation.}

\item{predError}{The predicted level of inaccuracy for the PSMs, obtained from external sources. Either the \code{ranger} model or \code{predError} must be specified.}

\item{ratioCalc}{A logical value indicating whether relative reporter intensities are calculated using the total reporter intensities across all channels.}
}
\value{
A data frame containing protein abundance estimates.
}
\description{
Aggregates PSMs using a random forest model.
}
\examples{
data(sample.PSM.FP)
regr <- loadPretrainedModel(useAvgCV=FALSE)
samples <- c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4',
    'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9')
groups <- samples
df <- getAttributes(
    PSM=sample.PSM.FP,
    fixedPTMs=c('229.1629', '57.0214'),
    colOfReporterIonInt=samples,
    groups=groups,
    setProgressBar=TRUE
)
aggregated_results <- aggregateByAttributes(
    PSM=df,
    colOfReporterIonInt=samples,
    ranger=regr,
    ratioCalc=FALSE
)
}
