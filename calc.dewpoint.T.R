calc.dewpoint.T = function(Substances=c(),
                           Fractions=c(),
                           Pressure=NULL,
                           UNIFAC=FALSE){
  # Load functions and define:
  source('sub.check.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('UNIFAC.gen.R');
  source('UNIFAC.R');
  source('calc.dewpoint.P.R');
  Ps = function(Temperature, ...) {
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature)
  } # list with saturated pressures
  if(UNIFAC){
    objF = function(Temperature, ...){
      abs(Pressure - calc.dewpoint.P(Substances,Fractions,Temperature)$Pressure);
    }
  }
  else{
    objF = function(Temperature, ...){ # objective function
      abs(sum(Fractions * Pressure / (Ps(Temperature))) - 1);
    }    
  }  
  e = 1e-7; # accuracy
  nos = length(Substances) # count number of substances 
  if(is.null(Pressure)){
  stop('No pressure [Pa] specified. Abort.');
  }
  if (nos == 0){
    stop('No substance specified. Abort.');
  }
  if (length(Fractions) != nos){ # check: # substances == # fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(sum(Fractions) - 1) > e){
    stop('Sum of mole fractions not equal to 1. Abort.');
  }
  if((!is.na(match(0, Fractions))) || (!is.na(match(-0, Fractions)))){
    stop('Negative or zero fraction found. Abort');    
  }
  sublist = NULL;
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i]);
    sublist = rbind(sublist, tmp[1,]);
    Tboil = Antoine.T(sublist$A[i],sublist$B[i],sublist$C[i],Pressure);
    sublist[i,7] = Tboil;
  }
  Tmin = min(sublist$T.boil);
  Tmax = max(sublist$T.boil);
  result = optimize(f=objF,lower=Tmin,upper=Tmax,tol=e);
  result = result[[1]]
  if(UNIFAC){
    tmp = calc.dewpoint.P(Substances,Fractions,result);
    result = list(Temperature=result,
                  ActivityCoefficients=tmp$Activity,
                  LiquidFractions=tmp$LiquidFractions);
  }
  return(result);
}