#' Get Average Scaled Error of log2FC for PSMs in a Spike-in Dataset
#' @description Calculates the Average Scaled Error of log2FC values required
#' for training sets.
#' @param PSM A data frame containing all PSMs used for training.
#' @param colOfReporterIonInt A vector of column names for reporter ion
#' intensities across different channels.
#' @param groups A vector specifying sample groups.
#' @param expectedRelativeAbundance A named list where group names are keys and
#' the corresponding expected relative abundance values for species at varying
#' concentrations are provided as values. Unknown ratios can be designated as
#' NA.
#' @param speciesAtConstLevel A string specifying the species that are spiked
#' in at a constant level.
#' @return A data frame containing PSMs with Average Scaled Error of log2FC
#' values required for the random forest model.
#' @importFrom stringr str_match
#' @importFrom utils combn
#' @examples
#' library(ExperimentHub)
#' library(stringr)
#' eh <- ExperimentHub()
#' benchmarkSet1 <- eh[['EH9637']]
#' # Load sample names (Sample 'H1+E1_1' ~ Sample 'H1+E6_3')
#' samples <- colnames(benchmarkSet1)[
#'     grep('H1[+]E[0-9]+_[1-4]', colnames(benchmarkSet1))
#' ]
#' groups <- str_match(samples, 'H1[+]E[0-9]+')[, 1]
#' PSM <- getAvgScaledErrorOfLog2FC(
#'     PSM=benchmarkSet1,
#'     colOfReporterIonInt=samples,
#'     groups=groups,
#'     expectedRelativeAbundance=list(`H1+E1`=1, `H1+E2`=2, `H1+E6`=NA),
#'     speciesAtConstLevel='HUMAN'
#' )
#' @export


getAvgScaledErrorOfLog2FC <- function(PSM, colOfReporterIonInt, groups,
                                      expectedRelativeAbundance,
                                      speciesAtConstLevel){

    PSM$Species <- str_match(PSM$Protein, '_([A-Z]+)$')[, 2]

    if(any(!is.na(expectedRelativeAbundance) & expectedRelativeAbundance < 0)){
        stop(
            'All expected relative abundance provided by ',
            'expectedRelativeAbundance should not be negative.'
        )
    }
    if(any(!groups %in% names(expectedRelativeAbundance))){
        stop('The expected relative abundance of some groups is not defined.')
    }

    expectedRelativeAbundance <- expectedRelativeAbundance[
        expectedRelativeAbundance != 0 & !is.na(expectedRelativeAbundance)
    ]

    cmps <- combn(names(expectedRelativeAbundance), m=2)
    abundance <- combn(unlist(expectedRelativeAbundance), m=2)

    for(idx in seq_len(ncol(cmps))){ # Get FC
        FC <- paste0('FC (', cmps[2, idx], ' vs ', cmps[1, idx], ')')
        PSM[, FC] <-
            apply(PSM[, colOfReporterIonInt[groups == cmps[2, idx]]], MARGIN =
                      1, FUN = mean) /
            apply(PSM[, colOfReporterIonInt[groups == cmps[1, idx]]], MARGIN =
                      1, FUN = mean)
    }

    for(idx in seq_len(ncol(cmps))){ # Get Error of log2FC
        expected_FC <- abundance[2, idx] / abundance[1, idx]
        FC <- paste0('FC (', cmps[2, idx], ' vs ', cmps[1, idx], ')')
        PSM[, paste0('Error of log2', FC)] <- abs(log2(PSM[, FC]) -
            ifelse(PSM$Species == speciesAtConstLevel, 0, log2(expected_FC)))
    }

    for(idx in seq_len(ncol(cmps))){ # Get Scaled Error of log2FC
        FC <- paste0('FC (', cmps[2, idx], ' vs ', cmps[1, idx], ')')
        for(j in unique(PSM$Species)){
            PSM[PSM$Species == j, paste0('Scaled Error of log2', FC)] <-
                PSM[PSM$Species == j, paste0('Error of log2', FC)] /
                mean(PSM[PSM$Species == j, paste0('Error of log2', FC)])
        }
    }

    # Get Average Scaled Error of log2FC
    if(length(grep(
        '^Scaled Error of log2FC [(].*?vs.*?[)]$',
        colnames(PSM)
    )) == 1){
        PSM$`Average Scaled Error of log2FC` <-
            PSM[, grep('^Scaled Error of log2FC [(].*?vs.*?[)]$',
                colnames(PSM))]
    }else{
        PSM$`Average Scaled Error of log2FC` <-
            apply(
                PSM[, grep('^Scaled Error of log2FC [(].*?vs.*?[)]$',
                    colnames(PSM))],
                MARGIN=1,
                FUN=mean
            )
    }
    return(PSM)
}
