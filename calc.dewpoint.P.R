calc.dewpoint.P = function(Substances=c(),
                           Fractions=c(),
                           Temperature=NULL){
  # function to calculate the dewpoint pressure of multicomponent mixtures
  # using UNIFAC to calculate the activity of the liquid phase
  # using no model to calculate the fugacity of the gas phase
  # load functions:
  source('sub.check.R');
  source('Antoine.P.R');
  source('UNIFAC.gen.R');
  source('UNIFAC.R');
  if(is.null(Temperature)){
    stop('No temperature [K] specified. Abort.');
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
  if (abs(sum(Fractions) - 1) > e){
    stop('Sum of mole fractions not equal to 1. Abort.');
  }
  if((!is.na(match(0, Fractions))) || (!is.na(match(-0, Fractions)))){
    stop('Negative or zero fraction found. Abort');    
  }
  sublist = NULL;
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i]);
    # maximum temperature check
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
  u = UNIFAC.gen(Substances); # load or generate UNIFAC values
  unu = u[[1]];
  aij = u[[2]];
  y = Fractions;
  phi = act = rep(1,nos);
  Ps =  Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature);
  c = 0;
  repeat{
    P = 1 / sum((phi * y) / (act * Ps));
    x = (phi * y * P) /(act * Ps);
    act = UNIFAC(x, unu, aij, Temperature);
    P = sum((act * x * Ps)/ phi);
    x = (phi * y * P) / (act * Ps);
    x = x / sum(x);
    act = UNIFAC(x, unu, aij, Temperature);
    Pnew = 1 / sum((phi * y) / (act * Ps));
    dP = abs(Pnew - P);
    c = c + 1;
    if((dP < e) || (c > 999))break
  }
  result = list(Pressure=Pnew,
                ActivityCoefficients=act,
                LiquidFractions=x,
                Iterations=c);
  return(result);
}