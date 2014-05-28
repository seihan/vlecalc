get.sub <- function(name){
  source('printf.R')
  subtbl <- read.csv('unifac-tbl3.csv',sep=',')
  sub <- subtbl[grep(paste('^',name,sep='',' '),subtbl$Name),]
  if (nrow(sub)==0){
    sub <- subtbl[grep(paste('^',name,sep='','$'),subtbl$Name),]
  }
  if (nrow(sub)==0){
    printf('\nSubstance "%s" not found.\n',name)
  }
  return(sub)
}