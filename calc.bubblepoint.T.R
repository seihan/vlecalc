calc.bubblepoint.T <- function(Substances=c(),
                               Fractions=c(),
                               Pressure=NULL,
                               UNIFAC=FALSE){
  # function to calculate the bubblepoint temperature of multicomponent mixtures
  # using UNIFAC to calculate the activity of the liquid phase
  # using no model to calculate the fugacity of the gas phase
  # load functions:
  source('Antoine.P.R')
  source('Antoine.T.R')
  source('sub.check.R')
  source('UNIFAC.gen.R')
  source('UNIFAC.R')
  Ps = function(Temperature, ...){
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature);
  }    
  checkTemp = function(Substances, Temperature, ...){
    # choose the right antoine coe
    sublist = NULL;
    for (i in 1:length(Substances)){
      tmp = sub.check(Substances[i]);
      c = 1;
      inrange = FALSE;
      repeat{
        if((Temperature > tmp$T.max[c] - 10) && ((c + 1) <= nrow(tmp))){
          c = c + 1;
        }
        else inrange = TRUE;
        if ((c == nrow(tmp)) || (inrange)) break;
      }
      if(Temperature > tmp$T.max[c]){
        warning('The temperature of ', round(Temperature,0),
                'K exceeds valid range of ', round(tmp$T.max[c],0), 'K from substance ', tmp[1,1],'.');
      }
      sublist = rbind(sublist, tmp[c,]);
    }
    return(sublist)
  }
  objF = function(Temperature, ...){
    if(!is.na(match(T,(Temperature > sublist$T.max)))){ # maximum temperature check
      sublist = checkTemp(Substances, Temperature);
    }
    Ps = Antoine.P(sublist$A,
                   sublist$B,
                   sublist$C,
                   Temperature);
    returnVal = abs(sum(Fractions * Ps)- Pressure);
    return(returnVal)
  }
  if(is.null(Pressure)){
    stop('No pressure [Pa] specified. Abort.');
  }
  nos = length(Substances); # count the number of substances
  e = 1e-7;  # accuracy
  if (nos == 0){
    stop('No substance specified. Abort.');
  }
  if (length(Fractions) != nos){ # check: # substances == # fractions
    stop('Number of substances differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(sum(Fractions) - 1) > e ){
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
  Tmin = min(sublist$T.boil); # lower limit
  Tmax = max(sublist$T.boil); # upper limit
  Temperature = optimize(f=objF,lower=Tmin,upper=Tmax,tol=e)[[1]];# ideal first guess
  result = Temperature;
  x = Fractions
  if(UNIFAC){
    u = UNIFAC.gen(Substances); # load or generate UNIFAC values
    unu = u[[1]];
    aij = u[[2]];
    kR = 8.314; # universal gas constant
    c = 0;  # iteration count
    x = Fractions;
    rat = 0;
    y = c();
    while (c < 1000) { # emergency break
      if(!is.na(match(T,(Temperature > sublist$T.max)))){ # maximum temperature check
        sublist = checkTemp(Substances, Temperature);
      }
      act = UNIFAC(x, unu, aij, Temperature);
      Pcalc = sum(act * Ps(Temperature) * x);  # calculate pressure
      rat = Pressure / Pcalc;
      if (abs(rat - 1) < e) break; # convergence criteria
      y =  Ps(Temperature) * x / Pcalc;
      v = kR * Temperature / Pcalc;
      f = act * Ps(Temperature) * exp((v * x) * (Pcalc - Pressure)/(kR * Temperature));
      PsNew = (y * act * Pressure / (x * f)) * Ps(Temperature);
      Temperature = sum(x * Antoine.T(sublist$A,
                                      sublist$B,
                                      sublist$C,
                                      PsNew));
      c = c + 1;  # count the iteration
    }
    y = y / sum(y);
    result = list(Temperature=Temperature,
                  ActivityCoefficients=act,
                  GasFractions=y,
                  Iterations=c);    
  }
  return(result);
}