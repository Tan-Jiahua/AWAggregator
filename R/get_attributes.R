#' Get attributes for PSMs
#' @description Retrieves attributes required for training or test sets.
#' @param PSM A data frame containing all PSMs used for training.
#' @param fixed_PTMs A numeric vector with the masses of fixed post-translational modifications (PTMs) in Da. Other PTMs will be treated as variable PTMs.
#' @param col_of_reporter_ion_int A vector of column names for reporter ion intensities across different channels.
#' @param groups A vector specifying sample groups.
#' @param groups_excluded_from_CV A vector of sample groups excluded from average CV calculations, which may occur due to a zero spike-in concentration for a species.
#' @param set_progress_bar A logical value indicating whether to display a progress bar.
#' @return A data frame containing the PSM table with attributes required for the random forest model.
#' @importFrom Peptides aaList
#' @importFrom stringr str_count
#' @importFrom stats median sd
#' @export


get_attributes = function(PSM, fixed_PTMs, col_of_reporter_ion_int, groups, groups_excluded_from_CV = NA, set_progress_bar = T){

  fixed_PTMs = paste0(fixed_PTMs, collapse = '|')
  PSM$Species = str_match(PSM$Protein, '_([A-Z,a-z, ]+)$')[, 2]

  # Get # variable modifications
  PSM$`# Variable Modifications` = str_count(PSM$`Assigned Modification`, '\\(.+?\\)') - str_count(PSM$`Assigned Modification`, fixed_PTMs)

  groups_for_CV = unique(groups[duplicated(groups) & !groups %in% groups_excluded_from_CV])

  if(!is.na(groups_excluded_from_CV)){
    print('These groups are removed when average CV is calculated because of the setting of groups_excluded_from_CV:')
    print(paste0(groups_excluded_from_CV, collapse = ', '))
  }
  if(any(!groups %in% groups_for_CV & !groups %in% groups_excluded_from_CV)){
    print('These groups are automatically removed when average CV is calculated because of lack of replicates:')
    print(paste0(unique(groups[!groups %in% groups_for_CV & !groups %in% groups_excluded_from_CV]), collapse = ', '))
  }

  # Scale hyperscore
  tmp = log2(PSM$Hyperscore)
  PSM$`Scaled Hyperscore` = (tmp - mean(tmp)) / sd(tmp)

  # Scale intensities
  tmp = ifelse(PSM$Intensity == 0, 0.5 * min(PSM[PSM$Intensity > 0, 'Intensity']), PSM$Intensity)
  tmp = log2(tmp)
  PSM$`Scaled Intensity` = (tmp - mean(tmp)) / sd(tmp)

  # Scale # variable modifications
  for(species in unique(PSM$Species)){
    PSM[PSM$Species == species, 'Scaled # Variable Modifications'] = PSM[PSM$Species == species, '# Variable Modifications'] - median(PSM[PSM$Species == species, '# Variable Modifications'])
  }

  if(length(groups_for_CV) != 0){
    # Get average CV
    tmp = rep(0, nrow(PSM))
    for(i in groups_for_CV){
      tmp = tmp + apply(PSM[, col_of_reporter_ion_int[groups == i]], MARGIN = 1, FUN = function(X) sd(X) / mean(X) * 100)
    }
    PSM$`Average CV [%]` = tmp / length(groups_for_CV)
  }else{
    print('There are no replicates so average CV will not be generated as an attribute.')
  }

  # Get distance metric
  PSM$`Distance Metric` = get_dist_metric(PSM, col_of_reporter_ion_int, set_progress_bar = set_progress_bar)
  PSM[is.na(PSM$`Distance Metric`), 'Distance Metric'] = mean(PSM$`Distance Metric`, na.rm = T)

  # Get AA composition fraction
  for(AA in aaList()){
    PSM[, paste0(AA, ' [%]')] = str_count(PSM[, 'Peptide'], AA) / nchar(PSM[, 'Peptide']) * 100
  }

  return(PSM)
}
