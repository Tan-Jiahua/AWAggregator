#' Benchmark Set 1 from FragPipe
#'
#' This data frame represents peptide spectrum matches (PSMs) in Benchmark Set 1, obtained from the search results of FragPipe. Columns unnecessary for the AWAggregator have been removed from the sample data.
#'
#' @format A data frame with 175293 rows and 21 variables:
#' \describe{
#'   \item{Peptide}{Peptide amino acid sequence}
#'   \item{Charge}{The charge state of the identified peptide}
#'   \item{Calibrated Observed Mass}{Mass of the identified peptide after m/z calibration (in Da)}
#'   \item{Calibrated Observed M/Z}{Mass-to-charge ratio of the peptide ion after m/z calibration}
#'   \item{Delta Mass}{Difference between calibrated observed peptide mass and theoretical peptide mass (in Da)}
#'   \item{Hyperscore}{The similarity score between observed and theoretical spectra}
#'   \item{Number of Missed Cleavages}{Number of potential enzymatic cleavage sites within the identified sequence}
#'   \item{Intensity}{Raw integrated precursor abundance for each PSM}
#'   \item{Assigned Modifications}{Post-translational modifications within the identified sequence}
#'   \item{Purity}{The proportion of total ion abundance in the inclusion window from the precursor}
#'   \item{Protein}{Protein sequence header corresponding to the identified peptide sequence}
#'   \item{H1+E1_1, H1+E1_2, H1+E1_3, H1+E2_1, H1+E2_2, H1+E2_3, H1+E2_4, H1+E6_1, H1+E6_2, H1+E6_3}{Processed reporter ion intensities for various samples. For example, "H1+E2_1" represents the reporter ion intensity for the first replicate (_1) of the group H1+E2, which contains a constant level of human proteins (H1) and a spiking of E. coli proteins (E2) at twice the base level}
#' }
#'
'benchmark.set.1'
