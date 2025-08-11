#' Sample PSM Data from Proteome Discoverer
#'
#' This data frame represents sample peptide spectrum matches (PSMs) mapped to
#' the proteins A0AV96, A0AVF1, A0AVT1, A0FGR8, and A0M8Q6, obtained from the
#' search results of Proteome Discoverer. Columns unnecessary for the
#' AWAggregator have been removed from the sample data.
#'
#' @docType data
#' @format A data frame with 128 rows and 21 variables:
#' \describe{
#'   \item{Annotated Sequence}{The names of the flanking residues of a peptide
#'   in a protein}
#'   \item{Modifications}{The static and dynamic modifications identified in
#'   the peptide}
#'   \item{Number of Proteins}{The number of mapped proteins}
#'   \item{Master Protein Accessions}{A description of the master proteins}
#'   \item{Number of Missed Cleavages}{The number of potential enzymatic
#'   cleavage sites within the identified sequence}
#'   \item{Charge}{The charge state of the peptide}
#'   \item{mz in Da}{The mass-to-charge ratio of the precursor ion, in daltons}
#'   \item{MHplus in Da}{The measured protonated monoisotopic mass of the
#'   peptides, in daltons}
#'   \item{Delta mz in Da}{The difference between the measured charged mass
#'   (m/z in Da) and the theoretical mass of the same charge (z)}
#'   \item{Isolation Interference in Percent}{The percentage of interference by
#'   co-isolation within the precursor isolation window}
#'   \item{Average Reporter SN}{The average reporter S/N values}
#'   \item{XCorr}{Scores the number of fragment ions that are common to two
#'   different peptides with the same precursor mass and calculates the
#'   cross-correlation score for all candidate peptides queried from the
#'   database}
#'   \item{Sample 1, Sample 2, Sample 3, Sample 4, Sample 5, Sample 6,
#'   Sample 7, Sample 8, Sample 9}{Processed reporter ion intensities from
#'   sample 1 to 9}
#' }
#' @usage data(sample.PSM.PD)
#'
'sample.PSM.PD'
