get.unifac <- function(formula){
  # function to get UNIFAC paramters from table
  # get.unifac('formula') -> Group, Symbol, Rk, Qk
  utbl <- read.csv('unifac-tbl1.csv',sep=',') # read UNIFAC surface paramters from csv file
  # search the UNIFAC table for functional group by formula and returns the line or 0
  # formula = case sensitive string
  line <- utbl[grep(paste('^',formula,sep='',' '),utbl$Symbol),]
  if (nrow(line)==0){
    line <- utbl[grep(paste('^',formula,sep='','$'),utbl$Symbol),]
  }
  if (nrow(line)==0){
    stop('\nGroup ',formula,' not found.\n')
  }
  return(line)
}