#' Aggregate PSMs to Protein Level
#' @description Aggregates PSMs using a random forest model.
#' @param PSM A data frame containing all PSMs to be aggregated.
#' @param col_of_reporter_ion_int A vector of column names representing reporter ion intensities across different channels.
#' @param ranger The random forest model to be applied for aggregation.
#' @param pred_error The predicted level of inaccuracy for the PSMs, obtained from external sources. Either the `ranger` model or `pred_error` must be specified.
#' @param ratio_calc A logical value indicating whether relative reporter intensities are calculated using the total reporter intensities across all channels.
#' @return A data frame containing protein abundance estimates.
#' @importFrom dplyr across group_by summarize %>%
#' @importFrom ranger ranger
#' @importFrom stats weighted.mean
#' @importFrom tidyr all_of
#' @export


aggregate_by_attributes = function(PSM, col_of_reporter_ion_int, ranger = NULL, pred_error = NULL, ratio_calc = F){

  if(missing(PSM)){
    stop('Variable "PSM" needs to be assigned')
  }
  if(missing(col_of_reporter_ion_int)){
    stop('Variable "col_of_reporter_ion_int" needs to be assigned')
  }
  if(is.null(ranger) & is.null(pred_error)){
    stop('Either variable "ranger" or variable "pred_error" needs to be assigned')
  }
  if(!is.null(ranger) & !is.null(pred_error)){
    stop('Variables "ranger" and "pred_error" should not be assigned simultaneously')
  }
  if(!is.null(pred_error) & length(pred_error) != nrow(PSM)){
    stop('The length of the vector "pred_error" should be the same as the number of rows of "PSM".')
  }
  if(!is.null(ranger)){
    PSMs_to_be_predicted = apply(PSM[, ranger$forest$independent.variable.names], MARGIN = 1, FUN = function(X) !any(is.na(X))) # Whether PSMs are predictable by random forest. These PSMs should not have any missing PSM attributes.
    if(any(!PSMs_to_be_predicted)){
      warning(paste0(sum(!PSMs_to_be_predicted), ' PSM(s) have (has) missing PSM attributes.'))
    }
    pred_error[PSMs_to_be_predicted] = predict(ranger, data = PSM[PSMs_to_be_predicted, ranger$forest$independent.variable.names])$predictions
    pred_error[PSMs_to_be_predicted & pred_error <= 0] = min(pred_error[PSMs_to_be_predicted & pred_error > 0]) # In case the prediction values are less than 0
    cnt1 = 0
    cnt2 = 0
    cnt3 = 0
    cnt4 = 0
    for(prot in PSM[!PSMs_to_be_predicted, 'Protein']){
      if(any(PSMs_to_be_predicted & PSM$Protein == prot)){
        pred_error[!PSMs_to_be_predicted & PSM$Protein == prot] = mean(pred_error[PSMs_to_be_predicted & PSM$Protein == prot])
        cnt1 = cnt1 + sum(!PSMs_to_be_predicted & PSM$Protein == prot)
        cnt2 = cnt2 + 1
      }else{
        pred_error[!PSMs_to_be_predicted & PSM$Protein == prot] = 1
        cnt3 = cnt3 + sum(!PSMs_to_be_predicted & PSM$Protein == prot)
        cnt4 = cnt4 + 1
      }
    }
    if(cnt1 != 0){
      warning(paste0(cnt2, ' protein(s) have (has) ', cnt1, ' PSM(s) with missing PSM attributes. The weight of these (this) PSM(s) is assigned as the average weights of the PSMs mapped to the same protein.'))
    }
    if(cnt3 != 0){
      warning(paste0(cnt4, ' protein(s) have (has) ', cnt3, ' PSM(s) but they (it) do(es) not have any PSMs without missing attributes. Their (its) abundance is estimated by the average of the corresponding PSM reporter ion intensities.'))
    }
  }
  if(ratio_calc){
    PSM[, col_of_reporter_ion_int] = t(apply(PSM[, col_of_reporter_ion_int], MARGIN = 1, FUN = function(X) X / sum(X)))
  }
  PSM[, col_of_reporter_ion_int] = log2(PSM[, col_of_reporter_ion_int])
  PSM$Weight = 1 / pred_error
  prot_norm = group_by(PSM, Protein) %>% summarize(across(all_of(col_of_reporter_ion_int), ~weighted.mean(., Weight, na.rm = T)))
  prot_norm[, col_of_reporter_ion_int] = 2^prot_norm[, col_of_reporter_ion_int]
  return(prot_norm)
}
