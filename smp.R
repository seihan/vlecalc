smp = function(Substances=c(),
               Bubblepoint=NULL,
               Dewpoint=NULL,
               Pressure=NULL){
  # Minimize the differantial between simulated and target boiling-, and
  # dew temperatures, by finding the optimal fractions of the entered
  # substances.
  #
  # Arguments:                                            
  #   Substances:     substances, c('ethanol', 'n-decane',...)
  #   Boiling Point:  target boiling temperature
  #   Dew Point:      target dew temperatur
  #   Pressure:       system pressure                                                        
  #
  # Returns:
  #   Vector with the molar fractions of the substances, (x1,x2,...,xn).
  #
  # Load and define functions:
  source('calc.bubblepoint.T.R')
  source('calc.dewpoint.T.R')
  source('Antoine.P.R')
  source('Antoine.T.R')
  #source('gamp.R')
  source('calc.boiling.T.R')
  source('sub.check.R')
  objF = function(x, ...){
    x = x / sum(x);
    if(abs(1 - sum(x)) > 1e-7) warning("\nThe sum of the fractions differ 1.");
    returnVal = abs(Tbc(s,x,P) - Tbr) + abs(Tdc(s,x,P) - Tdr);
    if(returnVal < 0.1){
      print(returnVal);
      print(x);
    }
    return(returnVal);
  }
  Tbc = function(s,x,P){calc.bubblepoint.T(s,x,P)};# Tb Bubblepoint
  Tdc = function(s,x,P){calc.dewpoint.T(s,x,P)};# Td Dewpoint
  Tbr = Bubblepoint;
  Tdr = Dewpoint;
  s = Substances;
  P = Pressure;
  nos = length(s);
  e = 1e-2; # accuracy
  sublist = NULL; # initiate list to sort substances by boiling point
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i]);
    sublist = rbind(sublist, tmp[1,]);
    sublist[i,7] = calc.boiling.T(Substances[i], P);
  }
  differ = max(sublist[,7]) - Tdr;
  if(differ < 0){
    stop('\nThe Dewpoint has a difference of ',round(differ,2),
         ' K to the boilingpoint of the highboiler.\nAbort blending.')
  }
  differ = Tbr - min(sublist[,7]);
  if(differ < 0){
    stop('\nThe Boilingpoint has a difference of ',round(differ,2),
         ' K to the boilingpoint of the lowboiler.\nAbort blending.')
  }
#   sublist <- sublist[order(sublist$T.boil),]
#   for (i in 1:nos){
#     Substances[i] <- strsplit(as.character(sublist$Substance[i]),' ')[[1]][1]
#   }
 
  ptm = proc.time();
  result = optim(rep(1,nos)/nos,objF);
  print(proc.time() - ptm);
  return(result);
}