#' Convert Proteome Discoverer Output to Required Input Format
#' @description Converts output from Proteome Discoverer into the input format required by AWAggregator.
#' @param PSM A data frame containing the PSM table from Proteome Discoverer to be converted.
#' @param protein A data frame containing the corresponding protein table from Proteome Discoverer.
#' @param colOfReporterIonInt A vector of column names for reporter ion intensities across different channels.
#' @return A data frame in the format required by AWAggregator.
#' @importFrom stringr str_match
#' @examples
#' data(sample.PSM.PD)
#' data(sample.prot.PD)
#' samples <- c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4',
#'     'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9')
#' df <- convertPDFormat(
#'     PSM=sample.PSM.PD,
#'     protein=sample.prot.PD,
#'     colOfReporterIonInt=samples
#' )
#' @export


convertPDFormat <- function(PSM, protein, colOfReporterIonInt){

    if(any(apply(PSM[, colOfReporterIonInt], MARGIN=1, FUN=function(X) all(is.na(X))))){
        warning('Some PSMs do not have any reporter ion intensities and they are removed.')
        PSM <- PSM[apply(PSM[, colOfReporterIonInt], MARGIN=1, FUN=function(X) any(!is.na(X))), ]
    }
    if(any(PSM$`Number of Proteins` != 1)){
        warning('Some PSMs are mapped to multiple proteins and they are removed.')
        PSM <- PSM[PSM$`Number of Proteins` == 1, ]
    }
    PSM$Peptide <- str_match(toupper(PSM$`Annotated Sequence`), '[.]([A-Z]+)')[, 2]
    colnames(PSM)[colnames(PSM) == 'Modifications'] <- 'Assigned Modifications'
    colnames(PSM)[colnames(PSM) == 'Master Protein Accessions'] <- 'Protein'
    PSM$Protein <- paste0(PSM$Protein, '_', str_match(protein[match(PSM$Protein, protein$Accession), 'Description'], 'OS=(.*?) [A-Z]+=')[, 2])
    colnames(PSM)[colnames(PSM) == 'mz in Da'] <- 'Calibrated Observed M/Z'
    PSM$`Calibrated Observed Mass` <- PSM$`MHplus in Da` - 1.008 # Deduct the mass of H+
    PSM$`Delta Mass` <- PSM$`Delta mz in Da` * PSM$Charge
    PSM$Purity <- 1 - PSM$`Isolation Interference in Percent` / 100
    PSM$Intensity <- PSM$`Average Reporter SN` # PSM attribute Intensity will go through a scaling in later steps
    if(any(is.na(PSM$Intensity))){
        warning('Missing intensity values (Average Reporter SN) are imputed as half of the minimal intensity.')
        PSM[is.na(PSM$Intensity), 'Intensity'] <- 0.5 * min(PSM[!is.na(PSM$Intensity), 'Intensity'])
    }
    PSM[is.na(PSM$Intensity)]
    colnames(PSM)[colnames(PSM) == 'XCorr'] <- 'Hyperscore'

    return(PSM)
}
