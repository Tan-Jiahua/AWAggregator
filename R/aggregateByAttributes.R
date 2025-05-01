#' Aggregate PSMs to Protein Level
#' @description Aggregates PSMs using a random forest model.
#' @param PSM A data frame containing all PSMs to be aggregated.
#' @param colOfReporterIonInt A vector of column names representing reporter ion intensities across different channels.
#' @param ranger The random forest model to be applied for aggregation.
#' @param predError The predicted level of inaccuracy for the PSMs, obtained from external sources. Either the `ranger` model or `predError` must be specified.
#' @param ratioCalc A logical value indicating whether relative reporter intensities are calculated using the total reporter intensities across all channels.
#' @return A data frame containing protein abundance estimates.
#' @importFrom dplyr across group_by summarize %>%
#' @importFrom ranger ranger
#' @importFrom rlang .data
#' @importFrom stats predict weighted.mean
#' @importFrom tidyr all_of
#' @examples
#' data(sample.PSM.FP)
#' regr <- loadPretrainedModel(useAvgCV=FALSE)
#' samples <- c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4',
#'     'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9')
#' groups <- samples
#' df <- getAttributes(
#'     PSM=sample.PSM.FP,
#'     fixedPTMs=c('229.1629', '57.0214'),
#'     colOfReporterIonInt=samples,
#'     groups=groups,
#'     setProgressBar=TRUE
#' )
#' aggregated_results <- aggregateByAttributes(
#'     PSM=df,
#'     colOfReporterIonInt=samples,
#'     ranger=regr,
#'     ratioCalc=FALSE
#' )
#' @export


aggregateByAttributes <- function(PSM, colOfReporterIonInt, ranger=NULL, predError=NULL, ratioCalc=FALSE){

    if (is.null(ranger) & is.null(predError)) {
        stop('Either variable "ranger" or variable "predError" needs to be assigned')
    }
    if (!is.null(ranger) & !is.null(predError)) {
        stop('Variables "ranger" and "predError" should not be assigned simultaneously')
    }
    if (!is.null(predError) & length(predError) != nrow(PSM)) {
        stop('The length of the vector "predError" should be the same as the number of rows of "PSM".')
    }
    if (!is.null(ranger)) {
        isPredictable <- apply(PSM[, ranger$forest$independent.variable.names], MARGIN=1, FUN=function(X) ! any(is.na(X))) # Whether PSMs are predictable by random forest. These PSMs should not have any missing PSM attributes.
        if(any(!isPredictable)){
            warning('Some PSMs have missing attributes.')
        }
        predError[isPredictable] <- predict(ranger, data=PSM[isPredictable, ranger$forest$independent.variable.names])$predictions
        predError[isPredictable & predError <= 0] <- min(predError[isPredictable & predError > 0]) # In case the prediction values are less than 0
        cnt1 <- 0
        cnt2 <- 0
        for(prot in PSM[!isPredictable, 'Protein']){
            if(any(isPredictable & PSM$Protein == prot)){
                predError[!isPredictable & PSM$Protein == prot] <- mean(predError[isPredictable & PSM$Protein == prot])
                cnt1 <- cnt1 + sum(!isPredictable & PSM$Protein == prot)
            }else{
                predError[!isPredictable & PSM$Protein == prot] <- 1
                cnt2 <- cnt2 + sum(!isPredictable & PSM$Protein == prot)
            }
        }
        if(cnt1 != 0){
            warning('Some proteins have PSMs with missing attributes. The weight of these PSMs is assigned as the average weights of the PSMs mapped to the same protein.')
        }
        if(cnt2 != 0){
            warning('Some proteins do not have any PSMs without missing attributes. Their abundance is estimated by the average of the corresponding PSM reporter ion intensities.')
        }
    }
    if(ratioCalc){
        PSM[, colOfReporterIonInt] <- t(apply(PSM[, colOfReporterIonInt], MARGIN=1, FUN=function(X) X / sum(X)))
    }
    PSM[, colOfReporterIonInt] <- log2(PSM[, colOfReporterIonInt])
    PSM$Weight <- 1 / predError
    res <- group_by(PSM, .data$Protein) %>% summarize(across(all_of(colOfReporterIonInt), ~weighted.mean(., Weight, na.rm=TRUE)))
    res[, colOfReporterIonInt] <- 2^res[, colOfReporterIonInt]
    return(res)
}
