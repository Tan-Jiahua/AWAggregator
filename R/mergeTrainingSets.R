#' Merge Training Sets
#' @description Extracts a similar number of PSMs from each input dataset and merges them into a single training set.
#' @param PSMList A named list where dataset names are keys and the corresponding data frames of PSMs used for training are values.
#' @param numPSMs The minimum number of PSMs to extract from each dataset for merging.
#' @param setProgressBar A logical value indicating whether to display a progress bar.
#' @return A data frame containing the merged PSM table from a subset of each input dataset.
#' @importFrom dplyr bind_rows mutate
#' @import progress
#' @importFrom purrr reduce
#' @importFrom tidyr gather
#' @importFrom toOrdinal toOrdinal
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' benchmarkSet1 <- eh[['EH9637']]
#' benchmarkSet2 <- eh[['EH9638']]
#' benchmarkSet3 <- eh[['EH9639']]
#' PSM <- mergeTrainingSets(
#'     PSMList=list(
#'         `Benchmark Set 1`=benchmarkSet1,
#'         `Benchmark Set 2`=benchmarkSet2,
#'         `Benchmark Set 3`=benchmarkSet3
#'     ),
#'     numPSMs=min(nrow(benchmarkSet1), nrow(benchmarkSet2), nrow(benchmarkSet3)),
#' )
#' @export


mergeTrainingSets <- function(PSMList, numPSMs, setProgressBar=TRUE){
    if(is.null(names(PSMList))){
        for(i in seq_len(length(PSMList))){
            message('The ', toOrdinal(i), ' PSM table in the list "PSMList" is not named and its name will be generated as "Dataset ', i, '".')
            names(PSMList)[i] <- paste0('Dataset ', i)
        }
    }
    for(i in seq_len(length(PSMList))){
        if(is.null(names(PSMList)[i]) | names(PSMList)[i] == ''){
            message('The ', toOrdinal(i), ' PSM table in the list "PSMList" is not named and its name will be generated as "Dataset ', i, '".')
            names(PSMList)[i] <- paste0('Dataset ', i)
        }
    }
    listOfSelectedPSMs <- vector(mode='list', length=length(PSMList))
    names(listOfSelectedPSMs) <- names(PSMList)

    if(setProgressBar){
        pb <- progress_bar$new(format='[:bar] :current/:total (:percent) elapsed :elapsed eta :eta', total=length(PSMList) * numPSMs)
    }
    for(i in seq_len(length(PSMList))){
        PSM <- PSMList[[i]]
        pickableProt <- unique(PSM$Protein)
        PSM <- PSM[PSM$Protein %in% pickableProt, ]
        selectedPSMs <- data.frame(matrix(ncol=ncol(PSM), nrow=0))
        colnames(selectedPSMs) <- colnames(PSM)

        while(nrow(selectedPSMs) < numPSMs){
            if(setProgressBar){
                pb$update((nrow(selectedPSMs) / numPSMs + i - 1) / length(PSMList))
            }
            pickedIdx <- sample(seq_len(length(pickableProt)), 1)
            pickedProt <- pickableProt[pickedIdx]
            pickableProt <- pickableProt[-pickedIdx]
            selectedPSMs <- rbind(selectedPSMs, PSM[PSM$Protein == pickedProt, ]) # Why rbind not bind_rows: https://stackoverflow.com/questions/42887217/difference-between-rbind-and-bind-rows-in-r
        }
        listOfSelectedPSMs[[i]] <- selectedPSMs
        if(setProgressBar){
            pb$update(i / length(PSMList))
        }
    }
    if(setProgressBar){
        pb$terminate()
    }

    col <- reduce(lapply(listOfSelectedPSMs, colnames), intersect)
    PSM <- bind_rows(
        lapply(names(listOfSelectedPSMs), FUN=function(X) mutate(listOfSelectedPSMs[[X]][, col], Dataset=X))
    )
    return(PSM)
}
