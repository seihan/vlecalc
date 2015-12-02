calc.bubblepoint.P = function(Substances,
                              Fractions,
                              Temperature){
  # Function to calculate the bubblepoint pressure
  # of multicomponent mixtures, using UNIFAC.
  # Author:  Hannes Seidel  Hannes.Seidel@b-tu.de
  # Load functions:
  source('Antoine.P.R');
  source('sub.check.R');
  source('UNIFAC.gen.R');
  source('UNIFAC.R');
  nos = length(Substances); # count number of substances
  if (length(Fractions) != nos){ # check: # substances == # fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistency (sum(xi = 1)?)
  if (abs(sum(Fractions) - 1) > 0.0001 ){
    stop('Sum of mole fractions not equal to 1. Abort.');
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
        warning('The temperature of ', Temperature,
                'K exceeds valid range of ', round(tmp$T.max[c],0), 'K from substance ', tmp[1,1],'.');
      }
      sublist = rbind(sublist, tmp[c,]);
    }
    return(sublist)
  }
  sublist = NULL
  sublist = checkTemp(Substances, Temperature);
  u = UNIFAC.gen(Substances); # load or generate UNIFAC values
  unu = u[[1]];
  aij = u[[2]];
  x = Fractions
  Ps = Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature);
  act = UNIFAC(x, unu, aij, Temperature);
  P = sum(act * x * Ps);
  y = act * x * Ps / P;
  result = list(pressure=P,
                activityCoefficients=act,
                vaporFractions=y);
  return(result);
}