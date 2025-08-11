#' Get attributes for PSMs
#' @description Retrieves attributes required for training or test sets.
#' @param PSM A data frame containing all PSMs used for training.
#' @param fixedPTMs A numeric vector with the masses of fixed
#' post-translational modifications (PTMs) in Da. Other PTMs will be treated as
#' variable PTMs.
#' @param colOfReporterIonInt A vector of column names for reporter ion
#' intensities across different channels.
#' @param groups A vector specifying sample groups.
#' @param groupsExcludedFromCV A vector of sample groups excluded from average
#' CV calculations, which may occur due to a zero spike-in concentration for a
#' species.
#' @param setProgressBar A logical value indicating whether to display a
#' progress bar.
#' @return A data frame containing the PSM table with attributes required for
#' the random forest model.
#' @importFrom Peptides aaList
#' @importFrom stringr str_count
#' @importFrom stats median sd
#' @examples
#' library(ExperimentHub)
#' library(stringr)
#' eh <- ExperimentHub()
#' benchmarkSet3 <- eh[['EH9639']]
#' # Load sample names (Sample 'H1+Y0_1' ~ Sample 'H1+Y10_2')
#' samples <- colnames(benchmarkSet3)[
#'     grep('H1[+]Y[0-9]+_[1-2]', colnames(benchmarkSet3))
#' ]
#' groups <- str_match(samples, 'H1[+]Y[0-9]+')[, 1]
#' PSM <- getPSMAttributes(
#'     PSM=benchmarkSet3,
#'     fixedPTM=c('304.2071', '125.0476'),
#'     colOfReporterIonInt=samples,
#'     groups=groups,
#'     groupsExcludedFromCV='H1+Y0'
#' )
#' @export


getPSMAttributes <- function(PSM, fixedPTMs, colOfReporterIonInt, groups,
                             groupsExcludedFromCV=NA, setProgressBar=TRUE){

    fixedPTMs <- paste0(fixedPTMs, collapse='|')
    PSM$Species <- str_match(PSM$Protein, '_([A-Z,a-z, ]+)$')[, 2]

    # Get # variable modifications
    PSM$`# Variable Modifications` <-
        str_count(PSM$`Assigned Modification`, '\\(.+?\\)') -
        str_count(PSM$`Assigned Modification`, fixedPTMs)

    groups_for_CV <- unique(
        groups[duplicated(groups) & !groups %in% groupsExcludedFromCV]
    )

    if(!is.na(groupsExcludedFromCV)){
        message(
            'These groups are removed when average CV is calculated because ',
            'of the setting of groupsExcludedFromCV:'
        )
        message(paste0(groupsExcludedFromCV, collapse=', '))
    }
    if(any(!groups %in% groups_for_CV & !groups %in% groupsExcludedFromCV)){
        message(
            'These groups are automatically removed when average CV is ',
            'calculated because of lack of replicates:'
        )
        message(paste0(unique(groups[
            !groups %in% groups_for_CV & !groups %in% groupsExcludedFromCV
        ]), collapse=', '))
    }

    # Scale hyperscore
    tmp <- log2(PSM$Hyperscore)
    PSM$`Scaled Hyperscore` <- (tmp - mean(tmp)) / sd(tmp)

    # Scale intensities
    tmp <- 0.5 * min(PSM[PSM$Intensity > 0, 'Intensity'])
    tmp <- ifelse(PSM$Intensity == 0, tmp, PSM$Intensity)
    tmp <- log2(tmp)
    PSM$`Scaled Intensity` <- (tmp - mean(tmp)) / sd(tmp)

    for(species in unique(PSM$Species)){ # Scale # variable modifications
        PSM[PSM$Species == species, 'Scaled # Variable Modifications'] <-
            PSM[PSM$Species == species, '# Variable Modifications'] -
            median(PSM[PSM$Species == species, '# Variable Modifications'])
    }

    # Get average CV
    if(length(groups_for_CV) != 0){
        tmp <- rep(0, nrow(PSM))
        for(i in groups_for_CV){
            tmp <- tmp +
                apply(
                    PSM[, colOfReporterIonInt[groups == i]],
                    MARGIN=1,
                    FUN=function(X) sd(X) / mean(X) * 100
                )
        }
        PSM$`Average CV [%]` <- tmp / length(groups_for_CV)
    }else{
        message(
            'There are no replicates so average CV will not be generated as ',
            'an attribute.'
        )
    }

    # Get distance metric
    PSM$`Distance Metric` <- getDistMetric(
        PSM,
        colOfReporterIonInt,
        setProgressBar=setProgressBar
    )
    PSM[is.na(PSM$`Distance Metric`), 'Distance Metric'] <-
        mean(PSM$`Distance Metric`, na.rm=TRUE)

    # Get AA composition fraction
    for(AA in aaList()){
        PSM[, paste0(AA, ' [%]')] <-
            str_count(PSM[, 'Peptide'], AA) / nchar(PSM[, 'Peptide']) * 100
    }

    return(PSM)
}
