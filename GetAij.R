source('printf.R')
get.aij <- function(i,j){
  atbl <- read.csv('unifac-tbl2.csv',sep=',') # read UNIFAC interaction paramters from csv file
  # search the UNIFAC interaction parameter table for the combination of groups
  for (k in 1:nrow(atbl)){
    if ((atbl[k,1] == i) && (atbl[k,2]== j)){
      return(atbl[k,3])
      break
    }
    else if (k == nrow(atbl)){
      printf('Combination of group %d and group %d not found.',i,j)
      return(0)
    }
  }
}