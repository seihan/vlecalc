smp = function(substances=NULL,
               fractions=NULL,
               pressure=NULL,
               samples=NULL,
               model=NULL){
  source('Substance.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('UNIFAC.R');
  source('UNIFAC.gen.R');
  source('SRK.R');
  Percent = function(x, digits = 0, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  PsAntoine = function(temperature, ...) {# vector with saturated pressures
    ps = Antoine.P(A, B, C, temperature);
    return(ps);
  }
  TsAntoine = function(pressure, ...){ # vector with saturated temperatures
    ts = Antoine.T(A, B, C, pressure);
    return(ts);
  }
  PsSRK = function(temperature, ...){
    ps = rep(0, nos);
    for(i in 1:nos){
      ps[i] = SRK(temperature=temperature,
                  x=1, Ac=Ac[i], Pc=Pc[i], Tc=Tc[i])$pressure;
    }
    return(ps)
  }
  TsSRK = function(pressure, ...){
    ts = rep(0, nos);
    for(i in 1:nos){
      ts[i] = SRK(pressure=pressure,
                  x=1, Ac=Ac[i], Pc=Pc[i], Tc=Tc[i])$temperature;
    }
    return(ts)
  }
  calc.bubble.T = function(x, ...){ # Antoine model
    ts = TsAntoine(pressure);
    temperature = sum(x * ts);
    c = 0; # iteration count
    while(c < 100){
      ps = PsAntoine(temperature);
      p = sum(x * ps);
      y = abs(x * ps / pressure); 
      rat = p / pressure;  # ratio of calculated and system pressure
      if(abs(1 - rat) < e){break;} # convergence criteria
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 100){
      printf('Warning: Iteration limit (Antoine). Calculated pressure = %3.3f Pa.\n', p); 
    } 
    result = list(temperature=temperature,
                  y=y,
                  pressure=p);
  }
  calc.bubbleUNIFAC.T = function(x, ...){ # Antoine UNIFAC
    ts = TsAntoine(pressure);
    temperature = sum(x * ts);
    c = 0;
    while(c < 100){
      act = UNIFAC(Fractions=x, unu=unu, aij=aij, Temperature=temperature); # activity coefficients
      ps = PsAntoine(temperature);
      p = sum(act * x * ps);
      y = act * x * ps / pressure; 
      rat = p / pressure; 
      if (abs(1 - rat) < e){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 100){
      printf('Warning: Iteration limit (Antoine). Calculated pressure = %3.3f Pa.\n', p); 
    }
    result = list(temperature=temperature,
                  y=abs(y))
    return(result);
  }
  calc.bubbleSRK.T = function(x, ...){ # cubic equation saturation pressure
    ts = TsSRK(pressure);
    temperature = sum(x * ts);
    c = 0;
    while(c < 100){
      ps = PsSRK(temperature);
      p = sum(x * ps);
      y = x * ps / p; 
      rat = p / pressure; 
      if (abs(1 - rat) < e){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 100) printf('Warning: Iteration limit (cubic, pure). Calculated pressure = %3.3f Pa.\n', pcalc)
    result = list(temperature=temperature,
                  y=abs(y));
    return(result);
  }
  calc.bubbleSRKU.T = function(x, ...){ # cubic saturation pressure + UNIFAC
    ts = TsSRK(pressure);
    temperature = sum(x * ts);
    c = 0;
    while(c < 100){
      act = UNIFAC(Fractions=x, unu=unu, aij=aij, Temperature=temperature);
      ps = PsSRK(temperature);
      p = sum(act * x * ps);
      y = act * x * ps / p; 
      rat = p / pressure; 
      if (abs(1 - rat) < e){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 100) printf('Warning: Iteration limit (cubic). Calculated pressure = %3.3f Pa.\n', pcalc)
    result = list(temperature=temperature,
                  y=abs(y));
    return(result);
  }
  printf = function(...) cat(sprintf(...));
  e = 1e-7;  # accuracy
  nos = length(substances); # count number of substances
  if(is.null(pressure)){
    stop('No pressure [Pa] specified. Abort.');
  }
  if(nos == 0){
    stop('No substance specified. Abort.');
  }
  if(is.null(fractions)){
    fractions = rep(1, nos) / nos;
  }
  if(length(fractions) != nos){ # check: N substances == N fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if(abs(1 - sum(fractions)) > e){
    stop('Sum of mole fractions not equal to 1. Abort.');
  }
  if(is.null(samples)) {
    stop('No samples to compare. Abort.')
  }
  if(is.null(model)){
    warning('No model declared, using model 1 of 5.')
    model = 1;
  } 
  Substances = rep(list(''), nos); # initiate the list of lists
  for(i in 1:nos){ # collect substances properties
    Substances[[i]] = Substance(substances[i]); # get properties
    A = Substances[[i]]$Antoine$A[1];
    B = Substances[[i]]$Antoine$B[1];
    C = Substances[[i]]$Antoine$C[1];
    tsat = Antoine.T(A, B, C, pressure); # calc boiling temperature using Antoine coefficients
    Substances[[i]]$Tsat = tsat;
    Ac = Substances[[i]]$Ac;
    Pc = Substances[[i]]$Pc;
    Tc = Substances[[i]]$Tc;
    TsatSRK = SRK(pressure, x=1, Ac=Ac, Pc=Pc, Tc=Tc); # calc boiling using cubic equation
    Substances[[i]]$TsatSRK = TsatSRK$temperature; 
  }
  names(Substances) = substances;
  Substances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  Suborg = Substances;
  abc=sapply(Substances, function(abc) abc$Antoine[1,1:3]);
  A = as.numeric(abc[1,]);
  B = as.numeric(abc[2,]);
  C = as.numeric(abc[3,]);
  Ac = sapply(Substances, function(ac) ac$Ac);
  Pc = sapply(Substances, function(pc) pc$Pc);
  Tc = sapply(Substances, function(tc) tc$Tc);
  substances = names(Substances)
  u = UNIFAC.gen(substances); # load UNIFAC values
  unu = u[[1]];
  aij = u[[2]];
  objF = function(fractions, 
                  plotting=F, ...){ # objective function for optimisation
    fractions = abs(fractions);
    fractions = fractions / sum(fractions);
    for(i in 1:nos){
      Substances[[i]]$MassFraction = fractions[i];
      Substances[[i]]$Mass = fractions[i] * 100;
      Substances[[i]]$Fraction = as.numeric(fractions[i] /
      Substances[[i]]$MolarMass); # temporary value
    }
    m = sum(sapply(Substances, function(mass) mass$Fraction)); 
    for(i in 1:nos){
      Substances[[i]]$MolarAmount = as.numeric(Substances[[i]]$Fraction / m);
      Substances[[i]]$Fraction = as.numeric(Substances[[i]]$MolarAmount);  
    }
    x = sapply(Substances, function(fractions) fractions$Fraction);
    mass = 100;
    temperatures = c();
    vapor = 0;
    while(round(mass, 0) > 0){
      if(model == 1){ # Antoine
        result = calc.bubble.T(x=x);
      }
      if(model == 2){ # Antoine UNIFAC
        result = calc.bubbleUNIFAC.T(x=x);
      }
      if(model == 3){ # cubic homogene
        result = SRK(pressure=pressure, x=x, Tc=Tc, Pc=Pc, Ac=Ac);
      } 
      if(model == 4){ # cubic sat pressure
        result = calc.bubbleSRK.T(x=x);
      }
      if(model == 5){ # cubic sat pressure + UNIFAC
        result = calc.bubbleSRKU.T(x=x);
      } 
      temperature = result$temperature;
      if(mass < 10){
        temp = length(temperatures);
        if(temperature < temperatures[temp]){ # check for raising temperature
          temperature = temperatures[temp];
        }
      }
      temperatures = c(temperatures, temperature);
      y = result$y;
      Mm = as.numeric(sapply(Substances, function(molmass) molmass$MolarMass));
      ym = (y * Mm) / sum(y * Mm); # convert molare to mass fraction
      for(i in 1:nos){
        Substances[[i]]$Mass = round(abs(as.numeric(Substances[[i]]$Mass - ym[i])), 9); # substract 1g vapor
      }
      masses = sapply(Substances, function(mass) mass$Mass);
      mass = sum(masses);
      xm = as.numeric(masses / mass);
      vapor = c(vapor, 100 - round(mass,0));
      x = (xm / Mm) / sum(xm / Mm);
      for(i in 1:nos){
        Substances[[i]]$Fraction = x[i];
      }
    }
    temp = length(temperatures)
    if((model == 1) || (model == 2)){
      temperatures = c(temperatures, Substances[[1]]$Tsat);      
    }
    else{
      temperatures = c(temperatures, Substances[[1]]$TsatSRK);      
    }
    boilingcurve = list(vapor=vapor,
                        temperature=temperatures);
    f = approxfun(boilingcurve$vapor, boilingcurve$temperature, method='linear');
    diff = sum(abs(f(samples$x) - samples$T)); # calc the difference between sample and simulation
    if(plotting){
      plot(boilingcurve$vapor, boilingcurve$temperature, col='green', type='l');
      points(samples$x, samples$T);
      return(boilingcurve)
    }
    print(diff)
    return(diff);
  }
  result = optim(par=fractions, fn=objF);
  boilingcurve = objF(fractions=result$par,
                      plotting=TRUE);
  result = list(substances=names(Substances),
                mass.fractions=abs(result$par)/sum(abs(result$par)),
                deviation=result$value,
                boilingcurve=boilingcurve);
  return(result);
}