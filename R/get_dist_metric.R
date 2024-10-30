#' Get Distance Metric
#' @description Calculates the distance metric for PSMs.
#' @param PSM A data frame containing the PSMs for which distance metrics are to be calculated.
#' @param channel A vector specifying the channels used for calculating the distance metric.
#' @param set_progress_bar A logical value indicating whether to display a progress bar.
#' @return A vector of distance metrics for the specified PSMs.
#' @import progress
#' @export

get_dist_metric = function(PSM, channel, set_progress_bar = T){
  PSM$info = paste0(PSM$Protein, ', ', PSM$Peptide)
  red = unique(PSM$info[duplicated(PSM$info)])
  PSM[PSM$info %in% red, 'group'] = PSM[PSM$info %in% red, 'info']
  red_PSM = unique(PSM$Protein[duplicated(PSM$Protein)])
  PSM[is.na(PSM$group) & PSM$Protein %in% red_PSM, 'group'] = PSM[is.na(PSM$group) & PSM$Protein %in% red_PSM, 'Protein']
  if(set_progress_bar){
    pb = progress_bar$new(format = '[:bar] :current/:total (:percent) elapsed :elapsed eta :eta', total = length(unique(PSM$group)))
  }
  mat = as.matrix(PSM[, channel])
  mat = t(apply(mat, MARGIN = 1, FUN = function(X) X / sum(X)))

  res = rep(NA, nrow(PSM))

  for(idx in unique(PSM$group)){
    if(set_progress_bar){
      pb$tick()
    }
    temp = which(PSM$group == idx)
    if(length(temp) == 1){# For single PSM in a group: mean distance to all other peptides of the protein
      res[temp] = mean(apply(mat[PSM$Protein == idx & PSM$group != idx, ], MARGIN = 1, FUN = function(X) sqrt(sum((mat[temp, ] - X)^2)))) # Edited
    }else if(length(temp) == 2){
      res[temp[1]] = res[temp[2]] = sqrt(sum((mat[temp[1], ] - mat[temp[2], ])^2))
    }else{
      res[temp] = sapply(temp, FUN = function(idx1){
        mean(sapply(temp[temp != idx1], FUN = function(idx2) sqrt(sum((mat[idx1, ] - mat[idx2, ])^2))))
      })
    }
  }
  res = as.numeric(res)
  return(res)
}
