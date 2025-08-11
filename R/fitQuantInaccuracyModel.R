#' Fit a Random Forest Model
#' @description Trains a random forest model to predict the level of
#' quantitative inaccuracy of PSMs.
#' @param PSM A data frame containing all PSMs used for training.
#' @param numTrees The number of trees to include in the random forest model.
#' @param useAvgCV A logical value indicating whether to include the average CV
#' attribute in the training. This parameter is ignored if `appliedAttributes`
#' is specified.
#' @param importance A logical value indicating whether to assess the
#' importance of attributes.
#' @param seed An integer seed for random number generation in the random
#' forest.
#' @param appliedAttributes A vector of attribute names to be used in training,
#' replacing the default attributes.
#' @return A trained random forest model.
#' @importFrom Peptides aaList
#' @importFrom ranger ranger
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
#'     # TMTpro tag (304.2071) and N-ethylmaleimide (125.0476) are applied as
#'     # fixed PTMs
#'     fixedPTM=c('304.2071', '125.0476'),
#'     colOfReporterIonInt=samples,
#'     groups=groups
#' )
#' PSM <- getAvgScaledErrorOfLog2FC(
#'     PSM=PSM,
#'     colOfReporterIonInt=samples,
#'     groups=groups,
#'     expectedRelativeAbundance=list(
#'         `H1+Y0`=0, `H1+Y1`=1, `H1+Y5`=5, `H1+Y10`=10
#'     ),
#'     speciesAtConstLevel='HUMAN'
#' )
#' regr <- fitQuantInaccuracyModel(PSM, useAvgCV=TRUE, seed=3979)
#' @export


fitQuantInaccuracyModel <- function(PSM, numTrees=500, useAvgCV=TRUE,
                                    importance=FALSE, seed,
                                    appliedAttributes=c(NA)){
    if(!is.na(appliedAttributes[1])){
        X <- PSM[, appliedAttributes]
    }else if(useAvgCV){
        X <- PSM[, c('Charge', 'Calibrated Observed Mass',
                     'Calibrated Observed M/Z', 'Delta Mass',
                     'Scaled Hyperscore', 'Number of Missed Cleavages',
                     'Scaled Intensity', 'Purity',
                     'Scaled # Variable Modifications', 'Average CV [%]',
                     'Distance Metric', paste0(aaList(), ' [%]'))]
    }else{
        X <- PSM[, c('Charge', 'Calibrated Observed Mass',
                     'Calibrated Observed M/Z', 'Delta Mass',
                     'Scaled Hyperscore', 'Number of Missed Cleavages',
                     'Scaled Intensity', 'Purity',
                     'Scaled # Variable Modifications', 'Distance Metric',
                     paste0(aaList(), ' [%]'))]
    }
    startTime <- Sys.time()
    regr <- ranger(num.trees=numTrees, seed=seed, importance='impurity', x=X,
                   y=PSM$`Average Scaled Error of log2FC`)
    endTime <- Sys.time()
    message(
        'Model training time = ',
        as.numeric(difftime(endTime, startTime, units = 'mins')),
        ' minutes'
    )
    return(regr)
}
