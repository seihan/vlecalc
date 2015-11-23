calc.dewpoint.T = function(Substances,
                           Fractions,
                           Pressure){
  # Load functions and define:
  source('sub.check.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('UNIFAC.gen.R');
  source('UNIFAC.R');
  source('calc.dewpoint.P.R');
  objF = function(Temperature, ...){
    abs(Pressure - calc.dewpoint.P(Substances,Fractions,Temperature)$Pressure);
  }
  e = 1e-7; # accuracy
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
  Temperature = result[[1]]
  tmp = calc.dewpoint.P(Substances,Fractions,Temperature);
  result = list(Temperature=Temperature,
                Activity=tmp$Activity,
                LiquidFractions=tmp$LiquidFractions);
  return(result);
}