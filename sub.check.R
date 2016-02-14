sub.check = function(name){ # check for exist of substance by name
#  subtbl = read.csv('Substances.Antoine-tbl.csv',sep=','); # load table
  substance = subtbl[grep(paste('^',name,sep='',' '),subtbl$Substance),];
  if (nrow(substance) == 0){
    substance = subtbl[grep(paste('^',name,sep='','$'),subtbl$Substance),];
  }
  if (nrow(substance) == 0){
    stop('\nSubstance "', name,'" not found.');
  }
  else{
    return(substance);
  }
}