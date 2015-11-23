calc.dewpoint.T = function(Substances, Fractions, Pressure){
  # function to calculate the dewpoint temperature of multicomponent mixtures
  # load functions:
  source('sub.check.R')
  source('Antoine.P.R')
  source('Antoine.T.R')
  source('calc.boiling.T.R')
  Ps = function(Temperature, ...) {
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature)
  } # list with saturated pressures
  kR = 8.314; # universal gas constant
  e = 1e-7;  # accuracy
  nos = length(Substances) # count number of substances
  if (length(Fractions) != nos){ # check: # substances == # fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(sum(Fractions) - 1) > e){
    stop('Sum of mole fractions not equal to 1. Abort.')
  }
  if((!is.na(match(0, Fractions))) || (!is.na(match(-0, Fractions)))){
    stop('Negative or zero fraction found. Abort');    
  }
  sublist = NULL;
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i]);
    sublist = rbind(sublist, tmp[1,]);
    sublist[i,7] = calc.boiling.T(Substances[i], Pressure);
  }
  y = Fractions;
  objF = function(Temperature, ...){ # objective function
    returnVal = abs(sum(y * Pressure / (Ps(Temperature))) - 1)
  }
  Tmin = min(sublist[,7]);
  Tmax = max(sublist[,7]);
  Temperature = optimize(f=objF,lower=Tmin,upper=Tmax,tol=e)[[1]];
  return(Temperature);
}