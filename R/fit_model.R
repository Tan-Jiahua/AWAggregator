#' Fit a Random Forest Model
#' @description Trains a random forest model to predict the level of quantitative inaccuracy of PSMs.
#' @param PSM A data frame containing all PSMs used for training.
#' @param num.trees The number of trees to include in the random forest model.
#' @param use_average_CV A logical value indicating whether to include the average CV attribute in the training. This parameter is ignored if `applied_attributes` is specified.
#' @param importance A logical value indicating whether to assess the importance of attributes.
#' @param seed An integer seed for random number generation in the random forest.
#' @param applied_attributes A vector of attribute names to be used in training, replacing the default attributes.
#' @return A trained random forest model.
#' @importFrom Peptides aaList
#' @importFrom ranger ranger
#' @export


fit_model = function(PSM, num.trees = 500, use_average_CV = T, importance = F, seed, applied_attributes = c(NA)){
  if(!is.na(applied_attributes[1])){
    X = PSM[, applied_attributes]
  }else if(use_average_CV){
    X = PSM[, c('Charge', 'Calibrated Observed Mass', 'Calibrated Observed M/Z', 'Delta Mass', 'Scaled Hyperscore', 'Number of Missed Cleavages', 'Scaled Intensity', 'Purity', 'Scaled # Variable Modifications', 'Average CV [%]', 'Distance Metric', paste0(aaList(), ' [%]'))]
  }else{
    X = PSM[, c('Charge', 'Calibrated Observed Mass', 'Calibrated Observed M/Z', 'Delta Mass', 'Scaled Hyperscore', 'Number of Missed Cleavages', 'Scaled Intensity', 'Purity', 'Scaled # Variable Modifications', 'Distance Metric', paste0(aaList(), ' [%]'))]
  }
  start_time = Sys.time()
  regr = ranger(num.trees = num.trees, seed = seed, importance = 'impurity', x = X, y = PSM$`Average Scaled Error of log2FC`)
  end_time = Sys.time()
  print(paste0('Model training time = ', as.numeric(difftime(end_time, start_time, units = 'mins')), ' minutes'))
  return(regr)
}
