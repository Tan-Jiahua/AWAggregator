#' Merge Training Sets
#' @description Extracts a similar number of PSMs from each input dataset and merges them into a single training set.
#' @param PSM_list A named list where dataset names are keys and the corresponding data frames of PSMs used for training are values.
#' @param num_PSMs The minimum number of PSMs to extract from each dataset for merging.
#' @param set_progress_bar A logical value indicating whether to display a progress bar.
#' @param seed An integer seed for random selection of PSMs.
#' @return A data frame containing the merged PSM table from a subset of each input dataset.
#' @importFrom dplyr bind_rows mutate
#' @import progress
#' @importFrom purrr reduce
#' @importFrom tidyr gather
#' @importFrom toOrdinal toOrdinal
#' @export


merge_training_sets = function(PSM_list, num_PSMs, set_progress_bar = T, seed){
  if(is.null(names(PSM_list))){
    for(i in 1:length(PSM_list)){
      warning(paste0('The ', toOrdinal(i), ' PSM table in the list "PSM_list" is not named and its name will be generated as "Dataset ', i, '".'))
      names(PSM_list)[i] = paste0('Dataset ', i)
    }
  }
  for(i in 1:length(PSM_list)){
    if(is.null(names(PSM_list)[i]) | names(PSM_list)[i] == ''){
      warning(paste0('The ', toOrdinal(i), ' PSM table in the list "PSM_list" is not named and its name will be generated as "Dataset ', i, '".'))
      names(PSM_list)[i] = paste0('Dataset ', i)
    }
  }
  list_of_selected_PSMs = vector(mode = 'list', length = length(PSM_list))
  names(list_of_selected_PSMs) = names(PSM_list)

  if(set_progress_bar){
    pb = progress_bar$new(format = '[:bar] :current/:total (:percent) elapsed :elapsed eta :eta', total = length(PSM_list) * num_PSMs)
  }
  for(i in 1:length(PSM_list)){
    PSM = PSM_list[[i]]
    pickable_prot = unique(PSM$Protein)
    PSM = PSM[PSM$Protein %in% pickable_prot, ]
    selected_PSMs = data.frame(matrix(ncol = ncol(PSM), nrow = 0))
    colnames(selected_PSMs) = colnames(PSM)

    set.seed(seed)
    while(nrow(selected_PSMs) < num_PSMs){
      if(set_progress_bar){
        pb$update((nrow(selected_PSMs) / num_PSMs + i - 1) / length(PSM_list))
      }
      picked_idx = sample(1:length(pickable_prot), 1)
      picked_prot = pickable_prot[picked_idx]
      pickable_prot = pickable_prot[-picked_idx]
      selected_PSMs = rbind(selected_PSMs, PSM[PSM$Protein == picked_prot, ]) # Why rbind not bind_rows: https://stackoverflow.com/questions/42887217/difference-between-rbind-and-bind-rows-in-r
    }
    list_of_selected_PSMs[[i]] = selected_PSMs
    if(set_progress_bar){
      pb$update(i / length(PSM_list))
    }
  }
  if(set_progress_bar){
    pb$terminate()
  }

  col = reduce(lapply(list_of_selected_PSMs, colnames), intersect)
  PSM = bind_rows(
    lapply(names(list_of_selected_PSMs), FUN = function(X) mutate(list_of_selected_PSMs[[X]][, col], Dataset = X))
  )
  return(PSM)
}
