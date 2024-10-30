#' Load a Pre-trained Random Forest Model
#' @description Loads a pre-trained random forest model for predicting the level of quantitative inaccuracy of PSMs.
#' @param use_average_CV A logical value indicating whether the average CV attribute was used in the training.
#' @return The pre-trained random forest model.
#' @export

load_pretrained_model = function(use_average_CV = T){
  if(use_average_CV){
    res = AWAggregator::regr
    res$forest$child.nodeIDs = AWAggregator::child.nodeIDs
    res$forest$split.varIDs = AWAggregator::split.varIDs
    res$forest$split.values = c(AWAggregator::split.values.1, AWAggregator::split.values.2)
  }else{
    res = AWAggregator::regr.no.CV
    res$forest$child.nodeIDs = AWAggregator::child.nodeIDs.no.CV
    res$forest$split.varIDs = AWAggregator::split.varIDs.no.CV
    res$forest$split.values = c(AWAggregator::split.values.no.CV.1, AWAggregator::split.values.no.CV.2)
  }
  return(res)
}
