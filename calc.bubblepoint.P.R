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
  sublist = NULL
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i]);
    sublist = rbind(sublist, tmp[1,]);
  }
  u = UNIFAC.gen(Substances); # load or generate UNIFAC values
  unu = u[[1]];
  aij = u[[2]];
  x = Fractions
  Ps = Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature);
  act = UNIFAC(x, unu, aij, Temperature);
  P = sum(act * x * Ps);
  y = act * x * Ps / P;
  return(P)
}