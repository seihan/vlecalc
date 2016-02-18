os.boilingline = function(substances=NULL,
                          fractions=NULL,
                          pressure=NULL){
  source('Substance.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('UNIFAC.R');
  source('UNIFAC.gen.R');
  source('SRK.R');
  Percent = function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  Ps = function(temperature, ...) {# list with saturated pressures
    A=sapply(Substances, function(Ant) Ant$Antoine$A[1]);
    B=sapply(Substances, function(Ant) Ant$Antoine$B[1]);
    C=sapply(Substances, function(Ant) Ant$Antoine$C[1]);
    fractions = sapply(Substances, function(frac) frac$Fraction);
    Ps = Antoine.P(A, B, C, temperature);
  }
  calc.bubble.T = function(...){
    while (c < 100) { # emergency break
      act = UNIFAC(x, unu, aij, temperature);
      ps = Ps(temperature)
      pcalc = sum(act * ps * x);  # calculate pressure
      rat = pressure / pcalc;
      if (abs(rat - 1) < e) break; # convergence criteria
      y =  ps * x / pcalc;
      v = kR * temperature / pcalc;
      f = act * ps * exp((v * x) * (pcalc - pressure)/(kR * temperature));
      PsNew = (y * act * pressure / (x * f)) * ps;
      A=sapply(Substances, function(Ant) Ant$Antoine$A[1]);
      B=sapply(Substances, function(Ant) Ant$Antoine$B[1]);
      C=sapply(Substances, function(Ant) Ant$Antoine$C[1]);
      temperature = sum(x * Antoine.T(A, B, C, PsNew));
      c = c + 1;  # count the iteration 
    }
    return(temperature);
  }
  calc.bubble.T.SRKUNIFAC = function(...){
    Ps = function(temperature, ...){
      ps = rep(0, nos);
      for(i in 1:nos){
        ps[i] = SRK(temperature=temperature, x=1, Ac=Ac[i], Pc=Pc[i], Tc=Tc[i])$pressure;
      }
      return(ps)
    }
    Ts = function(pressure, ...){
      ts = rep(0, nos);
      for(i in 1:nos){
        ts[i] = SRK(pressure=pressure[i], x=1, Ac=Ac[i], Pc=Pc[i], Tc=Tc[i])$temperature;
      }
      return(ts)
    }
    c = 0;
    while(c < 100){
      tempFormer = temperature;
      act = UNIFAC(x, unu, aij, temperature);
      ps = Ps(temperature);
      pcalc = sum(act * x * ps);
      y =  ps * x / pcalc;
      v = kR * temperature / pcalc;
      f = act * ps * exp((v * x) * (pcalc - pressure)/(kR * temperature));
      psNew = (y * act * pressure / (x * f)) * ps;
      temperature = sum(x * Ts(psNew));
      diff = abs(temperature - tempFormer);
      if (diff < e) break; # convergence criteria
      c = c + 1;  # count the iteration
    }
    return(temperature);
  }
  printf = function(...) cat(sprintf(...));
  ptm = proc.time(); # count the calc. time
  e = 1e-7;  # accuracy
  nos = length(substances); # count number of substances
  if(is.null(pressure)){
    stop('No pressure [Pa] specified. Abort.');
  }
  if (nos == 0){
    stop('No substance specified. Abort.');
  }
  if (length(fractions) != nos){ # check: N substances == N fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(1 - sum(fractions)) > e){
    stop('Sum of mole fractions not equal to 1. Abort.');
  }
  Substances = rep(list(''), nos); # initiate the list of lists
  for(i in 1:nos){
    Substances[[i]] = Substance(substances[i]); # get properties
    Substances[[i]]$Mass = fractions[i] * Substances[[i]]$MolarMass; # calc. the mass of Sub.
    Substances[[i]]$Fraction = fractions[i]; # add fractions
    A = Substances[[i]]$Antoine$A[1];
    B = Substances[[i]]$Antoine$B[1];
    C = Substances[[i]]$Antoine$C[1];
    Tsat = Antoine.T(A, B, C, pressure) # calc boiling temperature using Antoine coefficients
    Substances[[i]]$Tsat = Tsat;
    Ac = Substances[[i]]$Ac;
    Pc = Substances[[i]]$Pc;
    Tc = Substances[[i]]$Tc;
    TsatSRK = SRK(pressure, x=1, Ac=Ac, Pc=Pc, Tc=Tc); # calc boiling using cubic equation
    Substances[[i]]$TsatSRK = TsatSRK$temperature;
  }
  names(Substances) = substances;
  Substances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  print(Substances)
  m100 = sum(sapply(Substances, function(mass) mass$Mass)); # sum mass
  #printf('Mass = %3.3f g\n', m100);
  vapratio = 0;
  temperatures = c();
  temperaturesSRK = c();
  temperaturesSRKUNIFAC = c();
  while(nos > 1){
    temperature = Substances[[nos]]$Tsat;
    fractions = sapply(Substances, function(frac) frac$Fraction);
    Ac = sapply(Substances, function(ac) ac$Ac);
    Pc = sapply(Substances, function(pc) pc$Pc);
    Tc = sapply(Substances, function(tc) tc$Tc);
    substances = names(Substances)
    u = UNIFAC.gen(substances); # load UNIFAC values
    unu = u[[1]];
    aij = u[[2]];
    kR = 8.314; # universal gas constant
    c = 0;  # iteration count
    x = fractions;
    rat = 0;
    y = c();
    temperature = calc.bubble.T();
    temp.SRK = SRK(pressure=pressure, x=x, Tc=Tc, Pc=Pc, Ac=Ac)$temperature;
    temp.SRKUNIFAC = calc.bubble.T.SRKUNIFAC();
    printf('Iterations = %3.0f, Temperature = %3.2f, Temp.SRK = %3.2f, Temp.SRKUNIFAC = %3.2f\n',
           c, temperature, temp.SRK, temp.SRKUNIFAC);
    Substances[[nos]] = NULL; # remove the lowboiler
    nos = length(Substances);
    if(nos > 1){ # apply mass conservation
      m = sum(sapply(Substances, function(mass) mass$Mass)); # sum molar mass
      liquid = m / m100 * 100;
      n = rep(0, nos);
      for(i in 1:nos){
        xm = Substances[[i]]$Mass / m; # mass fraction
        n[i] = xm / Substances[[i]]$MolarMass; # molar amount
      }
      x = n / sum(n); # molar fraction
      for(i in 1:nos){
        Substances[[i]]$Fraction = x[i];
      }
    }
    else{
      liquid = 0;
    }
    vapor = 100 - liquid;
    vapratio = c(vapratio, vapor);
    temperatures = c(temperatures, temperature);
    temperaturesSRK = c(temperaturesSRK, temp.SRK);
    temperaturesSRKUNIFAC = c(temperaturesSRKUNIFAC, temp.SRKUNIFAC);
  }
  temperatures = c(temperatures, Substances[[1]]$Tsat);
  temperaturesSRK = c(temperaturesSRK, Substances[[1]]$TsatSRK);
  temperaturesSRKUNIFAC = c(temperaturesSRKUNIFAC, Substances[[1]]$TsatSRK);
  boilline = matrix(c(vapratio, temperatures), ncol=2)
  boillineSRK = matrix(c(vapratio, temperaturesSRK), ncol=2);
  boillineSRKUNIFAC = matrix(c(vapratio, temperaturesSRKUNIFAC), ncol=2);
  print(boilline)
  print(boillineSRK)
  print(boillineSRKUNIFAC)
  duration = proc.time() - ptm; # calculation time
  blp = polyfit(vapratio, temperatures, 3); # post processing -> polynomial regression
  blpSRK = polyfit(vapratio, temperaturesSRK, 3);
  blpSRKUNIFAC = polyfit(vapratio, temperaturesSRKUNIFAC, 3);
  #printf('\nCalculation takes %3.3f s', as.numeric(duration[1]))
  plot(vapratio, temperaturesSRK);
  points(vapratio, temperatures, pch=4);
  points(vapratio, temperaturesSRKUNIFAC, pch=2);
  lines(seq(0, 100, 1), polyval(p=blp, x=seq(0, 100, 1)), col='blue');
  lines(seq(0, 100, 1), polyval(p=blpSRK, x=seq(0, 100, 1)), col='green');
  lines(seq(0, 100, 1), polyval(p=blpSRKUNIFAC, x=seq(0, 100, 1)), col='red');  
  return(blp);
}