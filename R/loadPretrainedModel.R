#' Load a Pre-trained Random Forest Model
#' @description Loads a pre-trained random forest model for predicting the level of quantitative inaccuracy of PSMs.
#' @param useAvgCV A logical value indicating whether the average CV attribute was used in the training.
#' @return The pre-trained random forest model.
#' @importFrom ExperimentHub ExperimentHub
#' @examples
#' regr <- loadPretrainedModel(useAvgCV=TRUE)
#' @export

loadPretrainedModel <- function(useAvgCV=TRUE){
    eh <- ExperimentHub()
    if(useAvgCV){
        res <- eh[['EH9640']]
    }else{
        res <- eh[['EH9641']]
    }
    return(res)
}
