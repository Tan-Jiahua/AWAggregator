#' Get attributes for PSMs
#' @description Retrieves attributes required for training or test sets.
#' @param PSM A data frame containing all PSMs used for training.
#' @param fixedPTMs A numeric vector with the masses of fixed post-translational modifications (PTMs) in Da. Other PTMs will be treated as variable PTMs.
#' @param colOfReporterIonInt A vector of column names for reporter ion intensities across different channels.
#' @param groups A vector specifying sample groups.
#' @param groupsExcludedFromCV A vector of sample groups excluded from average CV calculations, which may occur due to a zero spike-in concentration for a species.
#' @param setProgressBar A logical value indicating whether to display a progress bar.
#' @return A data frame containing the PSM table with attributes required for the random forest model.
#' @importFrom Peptides aaList
#' @importFrom stringr str_count
#' @importFrom stats median sd
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' benchmarkSet3 <- eh[['EH9639']]
#' samples <- c('H1+Y0_1', 'H1+Y1_1', 'H1+Y5_1',
#'     'H1+Y10_1', 'H1+Y10_2', 'H1+Y5_2', 'H1+Y1_2', 'H1+Y0_2')
#' groups <- c('H1+Y0', 'H1+Y1', 'H1+Y5',
#'     'H1+Y10', 'H1+Y10', 'H1+Y5', 'H1+Y1', 'H1+Y0')
#' PSM <- getAttributes(
#'     PSM=benchmarkSet3,
#'     fixedPTM=c('304.2071', '125.0476'),
#'     colOfReporterIonInt=samples,
#'     groups=groups,
#'     groupsExcludedFromCV='H1+Y0'
#' )
#' @export


getAttributes <- function(PSM, fixedPTMs, colOfReporterIonInt, groups, groupsExcludedFromCV=NA, setProgressBar=TRUE){

    fixedPTMs <- paste0(fixedPTMs, collapse='|')
    PSM$Species <- str_match(PSM$Protein, '_([A-Z,a-z, ]+)$')[, 2]

    PSM$`# Variable Modifications` <- str_count(PSM$`Assigned Modification`, '\\(.+?\\)') - str_count(PSM$`Assigned Modification`, fixedPTMs) # Get # variable modifications

    groups_for_CV <- unique(groups[duplicated(groups) & !groups %in% groupsExcludedFromCV])

    if(!is.na(groupsExcludedFromCV)){
        message('These groups are removed when average CV is calculated because of the setting of groupsExcludedFromCV:')
        message(paste0(groupsExcludedFromCV, collapse=', '))
    }
    if(any(!groups %in% groups_for_CV & !groups %in% groupsExcludedFromCV)){
        message('These groups are automatically removed when average CV is calculated because of lack of replicates:')
        message(paste0(unique(groups[!groups %in% groups_for_CV & !groups %in% groupsExcludedFromCV]), collapse=', '))
    }

    tmp <- log2(PSM$Hyperscore) # Scale hyperscore
    PSM$`Scaled Hyperscore` <- (tmp - mean(tmp)) / sd(tmp)

    tmp <- ifelse(PSM$Intensity == 0, 0.5 * min(PSM[PSM$Intensity > 0, 'Intensity']), PSM$Intensity) # Scale intensities
    tmp <- log2(tmp)
    PSM$`Scaled Intensity` <- (tmp - mean(tmp)) / sd(tmp)

    for(species in unique(PSM$Species)){ # Scale # variable modifications
        PSM[PSM$Species == species, 'Scaled # Variable Modifications'] <- PSM[PSM$Species == species, '# Variable Modifications'] - median(PSM[PSM$Species == species, '# Variable Modifications'])
    }

    if(length(groups_for_CV) != 0){ # Get average CV
        tmp <- rep(0, nrow(PSM))
        for(i in groups_for_CV){
            tmp <- tmp + apply(PSM[, colOfReporterIonInt[groups == i]], MARGIN=1, FUN=function(X) sd(X) / mean(X) * 100)
        }
        PSM$`Average CV [%]` <- tmp / length(groups_for_CV)
    }else{
        message('There are no replicates so average CV will not be generated as an attribute.')
    }

    PSM$`Distance Metric` <- getDistMetric(PSM, colOfReporterIonInt, setProgressBar=setProgressBar) # Get distance metric
    PSM[is.na(PSM$`Distance Metric`), 'Distance Metric'] <- mean(PSM$`Distance Metric`, na.rm=TRUE)

    for(AA in aaList()){ # Get AA composition fraction
        PSM[, paste0(AA, ' [%]')] <- str_count(PSM[, 'Peptide'], AA) / nchar(PSM[, 'Peptide']) * 100
    }

    return(PSM)
}
