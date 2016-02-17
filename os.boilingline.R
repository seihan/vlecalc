os.boilingline = function(substances=NULL,
                          fractions=NULL,
                          pressure=NULL){
  source('Substance.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('SRK.R');
  Percent = function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  Ps = function(temperature, ...) {# list with saturated pressures
    A=sapply(Substances, function(Ant) Ant$Antoine$A[1]);
    B=sapply(Substances, function(Ant) Ant$Antoine$B[1]);
    C=sapply(Substances, function(Ant) Ant$Antoine$C[1]);
    fractions = sapply(Substances, function(frac) frac$fraction);
    Ps = Antoine.P(A, B, C, temperature);
  }
  printf = function(...) cat(sprintf(...));
  objF = function(temperature, ...){
    returnVal = abs(sum(fractions * Ps(temperature))- pressure);
    return(returnVal)
  }
  ptm = proc.time();
  e <- 1e-7  # accuracy
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
    Tsat = Antoine.T(A, B, C, pressure) # calc boiling temperature
    Substances[[i]]$Tsat = Tsat;
  }
  names(Substances) = substances;
  Substances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  print(Substances)
  m100 = sum(sapply(Substances, function(mass) mass$Mass)); # sum mass
  printf('Mass = %3.3f g\n', m100);
  vapratio = 0;
  temperatures = c();
  while(nos > 1){
    Tmin = Substances[[nos]]$Tsat; # lower limit
    #Tmax = Substances[[1]]$Tsat; # upper limit
    # temperature = optimize(f=objF, lower=Tmin, upper=Tmax, tol=e)[[1]];# ideal first guess
    temperature = Tmin
    fractions = sapply(Substances, function(frac) frac$Fraction);
    #fractions = fractions / sum(fractions)
    substances = names(Substances)
    u = UNIFAC.gen(substances); # load or generate UNIFAC values
    unu = u[[1]];
    aij = u[[2]];
    kR = 8.314; # universal gas constant
    c = 0;  # iteration count
    x = fractions;
    rat = 0;
    y = c();
    while (c < 100) { # emergency break
      act = UNIFAC(x, unu, aij, temperature);
      pcalc = sum(act * Ps(temperature) * x);  # calculate pressure
      rat = pressure / pcalc;
      if (abs(rat - 1) < e) break; # convergence criteria
      y =  Ps(temperature) * x / pcalc;
      v = kR * temperature / pcalc;
      f = act * Ps(temperature) * exp((v * x) * (pcalc - pressure)/(kR * temperature));
      PsNew = (y * act * pressure / (x * f)) * Ps(temperature);
      A=sapply(Substances, function(Ant) Ant$Antoine$A[1]);
      B=sapply(Substances, function(Ant) Ant$Antoine$B[1]);
      C=sapply(Substances, function(Ant) Ant$Antoine$C[1]);
      temperature = sum(x * Antoine.T(A, B, C, PsNew));
      c = c + 1;  # count the iteration 
    }
    printf('Iterations = %3.0f, Temperature = %3.2f\n', c, temperature);
    Substances[[nos]] = NULL; # remove the lowboiler
    nos = length(Substances);
    if(nos > 1){
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
    print(Substances)
  }
  temperatures = c(temperatures, Substances[[1]]$Tsat);
  boilline = matrix(c(vapratio, temperatures), ncol=2)
  duration = proc.time() - ptm; # calculation time
  print(boilline)
  plot(vapratio, temperatures)
  p = polyfit(vapratio, temperatures, 3); # post processing -> polynomial regression
  xf = seq(0, 100, 1);
  yf = polyval(p, xf);
  lines(xf, yf, col="red");
  print(duration[1])
}