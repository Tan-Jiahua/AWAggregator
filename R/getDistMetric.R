#' Get Distance Metric
#' @description Calculates the distance metric for PSMs. Distance metric
#' reflects on whether the quantified ratio of each pair of samples of a PSM
#' diverges from other PSMs in the same redundant/unique group. Redundant
#' group, unique group and distance metric were originally defined in the iPQF
#' method. Please refer to "iPQF: a new peptide-to-protein summarization method
#' using peptide spectra characteristics to improve protein quantification" for
#' more details.
#' @param PSM A data frame containing the PSMs for which distance metrics are
#' to be calculated.
#' @param channel A vector specifying the channels used for calculating the
#' distance metric.
#' @param setProgressBar A logical value indicating whether to display a
#' progress bar.
#' @return A vector of distance metrics for the specified PSMs.
#' @import progress
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' benchmarkSet3 <- eh[['EH9639']]
#' # Load sample names (Sample 'H1+Y0_1' ~ Sample 'H1+Y10_2')
#' samples <- colnames(benchmarkSet3)[
#'     grep('H1[+]Y[0-9]+_[1-2]', colnames(benchmarkSet3))
#' ]
#' df <- getDistMetric(
#'     PSM=benchmarkSet3,
#'     channel=samples,
#'     setProgressBar=TRUE
#' )
#' @export

getDistMetric <- function(PSM, channel, setProgressBar=TRUE){
    PSM$info <- paste0(PSM$Protein, ', ', PSM$Peptide)
    red <- unique(PSM$info[duplicated(PSM$info)])
    PSM[PSM$info %in% red, 'group'] <- PSM[PSM$info %in% red, 'info']
    red_PSM <- unique(PSM$Protein[duplicated(PSM$Protein)])
    PSM[is.na(PSM$group) & PSM$Protein %in% red_PSM, 'group'] <-
        PSM[is.na(PSM$group) & PSM$Protein %in% red_PSM, 'Protein']
    if(setProgressBar){
        pb <- progress_bar$new(
            format=
                '[:bar] :current/:total (:percent) elapsed :elapsed eta :eta',
            total=length(unique(PSM$group))
        )
    }
    mat <- as.matrix(PSM[, channel])
    mat <- t(apply(mat, MARGIN=1, FUN=function(X) X / sum(X)))

    res <- rep(NA, nrow(PSM))

    for(idx in unique(PSM$group)){
        if(setProgressBar){
            pb$tick()
        }
        tmp <- which(PSM$group == idx)
        # For single PSM in a group:
        # mean distance to all other peptides of the protein
        if(length(tmp) == 1){
            res[tmp] <- mean(apply(
                mat[PSM$Protein == idx & PSM$group != idx, ],
                MARGIN=1,
                FUN=function(X) sqrt(sum((mat[tmp, ] - X)^2))
        ))
        }else if(length(tmp) == 2){
            res[tmp[1]] <-
                res[tmp[2]] <- sqrt(sum((mat[tmp[1], ] - mat[tmp[2], ])^2))
        }else{
            res[tmp] <- vapply(tmp, FUN=function(idx1){
                mean(vapply(
                    tmp[tmp != idx1],
                    FUN=function(idx2) sqrt(sum((mat[idx1, ] - mat[idx2, ])^2)),
                    FUN.VALUE=numeric(1))
                )
            }, FUN.VALUE=numeric(1))
        }
    }
    res <- as.numeric(res)
    return(res)
}
