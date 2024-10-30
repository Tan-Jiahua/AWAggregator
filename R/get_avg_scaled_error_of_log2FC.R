#' Get Average Scaled Error of log2FC for PSMs in a Spike-in Dataset
#' @description Calculates the Average Scaled Error of log2FC values required for training sets.
#' @param PSM A data frame containing all PSMs used for training.
#' @param col_of_reporter_ion_int A vector of column names for reporter ion intensities across different channels.
#' @param groups A vector specifying sample groups.
#' @param expected_relative_abundance A named list where group names are keys and the corresponding expected relative abundance values for species at varying concentrations are provided as values. Unknown ratios can be designated as NA.
#' @param species_at_const_level A string specifying the species that are spiked in at a constant level.
#' @return A data frame containing PSMs with Average Scaled Error of log2FC values required for the random forest model.
#' @importFrom stringr str_match
#' @importFrom utils combn
#' @export


get_avg_scaled_error_of_log2FC = function(PSM, col_of_reporter_ion_int, groups, expected_relative_abundance, species_at_const_level){

  PSM$Species = str_match(PSM$Protein, '_([A-Z]+)$')[, 2]

  if(any(!is.na(expected_relative_abundance) & expected_relative_abundance < 0)){
    stop('All expected relative abundance provided by expected_relative_abundance should not be negative.')
  }
  if(any(!groups %in% names(expected_relative_abundance))){
    stop('The expected relative abundance of some groups is not defined.')
  }

  expected_relative_abundance = expected_relative_abundance[expected_relative_abundance != 0 & !is.na(expected_relative_abundance)]

  cmps = combn(names(expected_relative_abundance), m = 2)
  abundance = combn(unlist(expected_relative_abundance), m = 2)

  # Get FC
  for(idx in 1:ncol(cmps)){
    FC = paste0('FC (', cmps[2, idx], ' vs ', cmps[1, idx], ')')
    PSM[, FC] =
      apply(PSM[, col_of_reporter_ion_int[groups == cmps[2, idx]]], MARGIN = 1, FUN = mean) /
      apply(PSM[, col_of_reporter_ion_int[groups == cmps[1, idx]]], MARGIN = 1, FUN = mean)
  }

  # Get Error of log2FC
  for(idx in 1:ncol(cmps)){
    expected_FC = abundance[2, idx] / abundance[1, idx]
    FC = paste0('FC (', cmps[2, idx], ' vs ', cmps[1, idx], ')')
    PSM[, paste0('Error of log2', FC)] = abs(log2(PSM[, FC]) - ifelse(PSM$Species == species_at_const_level, 0, log2(expected_FC)))
  }

  # Get Scaled Error of log2FC
  for(idx in 1:ncol(cmps)){
    FC = paste0('FC (', cmps[2, idx], ' vs ', cmps[1, idx], ')')
    for(j in unique(PSM$Species)){
      PSM[PSM$Species == j, paste0('Scaled Error of log2', FC)] = PSM[PSM$Species == j, paste0('Error of log2', FC)] / mean(PSM[PSM$Species == j, paste0('Error of log2', FC)])
    }
  }

  # Get Average Scaled Error of log2FC
  if(length(grep('^Scaled Error of log2FC [(].*?vs.*?[)]$', colnames(PSM))) == 1){
    PSM$`Average Scaled Error of log2FC` = PSM[, grep('^Scaled Error of log2FC [(].*?vs.*?[)]$', colnames(PSM))]
  }else{
    PSM$`Average Scaled Error of log2FC` = apply(PSM[, grep('^Scaled Error of log2FC [(].*?vs.*?[)]$', colnames(PSM))], MARGIN = 1, FUN = mean)
  }
  return(PSM)
}
