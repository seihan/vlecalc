cs.boilingline = function(substances=NULL,
                          fractions=NULL,
                          pressure=NULL,
                          verbose=F){
  source('Substance.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('UNIFAC.R');
  source('UNIFAC.gen.R');
  source('SRK.R');
  Percent = function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  Antoine.limit = function(temperature, verbose=FALSE, ...){
    A = rep(0, nos);
    B = rep(0, nos);
    C = rep(0, nos);
    for(i in 1:nos){
      Arows = nrow(Substances[[i]]$Antoine);
      for(j in 1:Arows){
        tmax = Substances[[i]]$Antoine$T.max[j];
        if(temperature < tmax){
          A[i] = Substances[[i]]$Antoine$A[j];
          B[i] = Substances[[i]]$Antoine$B[j];
          C[i] = Substances[[i]]$Antoine$C[j];
          break;
        }
        else if((j == Arows) && (temperature > tmax)){
          A[i] = Substances[[i]]$Antoine$A[Arows];
          B[i] = Substances[[i]]$Antoine$B[Arows];
          C[i] = Substances[[i]]$Antoine$C[Arows];
          if(verbose){
            warning('Temperature ', round(temperature, 2),
                    ' K exceed validity range ',
                    Substances[[i]]$Antoine$T.max[j],
                    ' K of Substance ',
                    names(Substances[i]));
          }
        }
      }
    }
    return(list(A=A,
                B=B,
                C=C));
  }
  PsAntoine = function(temperature, ...) {# list with saturated pressures
    tmax = min(sapply(Substances, function(tmax) tmax$Antoine$T.max[1]));
    if(tmax < temperature){
      ABC = Antoine.limit(temperature);
      A = ABC$A;
      B = ABC$B;
      C = ABC$C;
    }
    else{
      A=sapply(Substances, function(Ant) Ant$Antoine$A[1]);
      B=sapply(Substances, function(Ant) Ant$Antoine$B[1]);
      C=sapply(Substances, function(Ant) Ant$Antoine$C[1]);
    }
    fractions = sapply(Substances, function(frac) frac$Fraction);
    ps = Antoine.P(A, B, C, temperature);
  }
  TsAntoine = function(pressure, ...){
    tmax = min(sapply(Substances, function(tmax) tmax$Antoine$T.max[1]));
    if(tmax < temperature){
      ABC = Antoine.limit(temperature);
      A = ABC$A;
      B = ABC$B;
      C = ABC$C;
    }
    else{
      A=sapply(Substances, function(Ant) Ant$Antoine$A[1]);
      B=sapply(Substances, function(Ant) Ant$Antoine$B[1]);
      C=sapply(Substances, function(Ant) Ant$Antoine$C[1]); 
    }
    ts = Antoine.T(A, B, C, pressure);
  }
  SRK.limit = function(temperature, ...){
    tmax = min(sapply(Substances, function(tmax) tmax$Tc))
    if(tmax < temperature){
      for(i in 1:nos){
        if(Substances[[i]]$Tc < temperature){
          warning('Temperature ', round(temperature, 2),
                  ' K exceed the critical temperature ',
                  Substances[[i]]$Tc,
                  ' K of Substance ',
                  names(Substances[i]));
        }
      }
    }
  }
  PsSRK = function(temperature, ...){
    ps = rep(0, nos);
    for(i in 1:nos){
      ps[i] = SRK(temperature=temperature, method = 'ps',
                  x=1, Ac=Ac[i], Pc=Pc[i], Tc=Tc[i])$pressure;
    }
    return(ps)
  }
  TsSRK = function(pressure, ...){
    ts = rep(0, nos);
    for(i in 1:nos){
      ts[i] = SRK(pressure=pressure, method = 'ts',
                  x=1, Ac=Ac[i], Pc=Pc[i], Tc=Tc[i])$temperature;
    }
    return(ts)
  }
  calc.bubble.T = function(...){
    ts = sapply(Substances, function(tsat) tsat$Tsat);
    temperature = sum(x * ts);
    c = 0;
    repeat{
      ps = PsAntoine(temperature);
      p = sum(x * ps);
      y = x * ps / p; 
      rat = p / pressure;
      if ((abs(rat - 1) < e) || (c > 99)){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 100){
      printf('Warning: Iteration limit (Antoine). Calculated pressure = %3.3f Pa.\n', p); 
    }
    return(temperature);
  }
  calc.bubbleSRK.T = function(...){
    ts = TsSRK(pressure);
    temperature = sum(x * ts);
    c = 0;
    repeat{
      ps = PsSRK(temperature);
      p = sum(x * ps);
      y = x * ps / p; 
      rat = p / pressure; 
      if ((abs(rat - 1) < e) || (c > 99)){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 100) printf('Warning: Iteration limit (cubic). Calculated pressure = %3.3f Pa.\n', pcalc)
    return(temperature);
  }
  calc.dewpoint.T = function(...){
    temperature = Substances[[1]]$Tsat;
    c = 0;
    repeat{
      K = PsAntoine(temperature) / pressure; # Equilibrium constants
      x = fractions / K;
      S = sum(x);
      if ((abs(S - 1) < 1e-4) || (c > 99)){break;}
      else{
        c = c + 1;
        temp = 0.01 * temperature * (S - 1) / S;
        temperature = temperature + temp;
        x = x / S;
        x = x / sum(x);
      }
    }
    return(list(temperature=temperature,
                x=x));
  } # end calc.dewpoint.T
  calc.dewpointSRK.T = function(...){
    temperature = Substances[[1]]$TsatSRK;
    c = 0;
    repeat{
      K = PsSRK(temperature) / pressure; # Equilibrium constants
      x = fractions / K;
      S = sum(x);
      if ((abs(S - 1) < 1e-4) || (c > 99)){break;}
      else{
        c = c + 1;
        temp = 0.01 * temperature * (S - 1) / S;
        temperature = temperature + temp;
        x = x / S;
        x = x / sum(x);
      }
    }
    return(list(temperature=temperature,
                x=x));
  }
  objF = function(l, X=FALSE, ...){ # function for optimization
    xf = l * (1 - K) + K;  # x factors for linear equation system
    A = matrix(0, nos, nos);  # initiate matrix
    for (i in 1:nos) {
      A[i, i] = xf[i]; # fill the matrix
    }
    x = solve(A, b);  # solve the system
    if(X){
      return(x);
    }
    else{
      return(abs(sum(x) - 1) + abs(sum(K * x) - 1));
    }
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
  for(i in 1:nos){ # collect substances properties
    Substances[[i]] = Substance(substances[i]); # get properties
    Substances[[i]]$MassFraction = fractions[i];
    Substances[[i]]$Fraction = fractions[i] / Substances[[i]]$MolarMass; # temporary value
    Arows = nrow(Substances[[i]]$Antoine);
    for(j in 1:Arows){
      A = Substances[[i]]$Antoine$A[j];
      B = Substances[[i]]$Antoine$B[j];
      C = Substances[[i]]$Antoine$C[j];
      tsat = Antoine.T(A, B, C, pressure); # calc boiling temperature using Antoine coefficients
      tmax = Substances[[i]]$Antoine$T.max[j];
      if(tsat < tmax){
        break;
      }
      if((tsat > tmax) && (j == Arows)){
        warning('Saturation temperature ', round(tsat, 2),
                ' exceed the validity of Antoine coefficients from substance ', substances[i], '.')
      }
    }
    Substances[[i]]$Tsat = tsat;
    Ac = Substances[[i]]$Ac;
    Pc = Substances[[i]]$Pc;
    Tc = Substances[[i]]$Tc;
    TsatSRK = SRK(pressure, x=1, Ac=Ac, Pc=Pc, Tc=Tc, method = 'ts'); # calc boiling using cubic equation
    Substances[[i]]$TsatSRK = TsatSRK$temperature; 
  }
  m = sum(sapply(Substances, function(mass) mass$Fraction)); 
  for(i in 1:nos){
    Substances[[i]]$MolarAmount = Substances[[i]]$Fraction / m;
    Substances[[i]]$Fraction = Substances[[i]]$MolarAmount;
  }
  names(Substances) = substances;
  Substances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  print(Substances)
  temperatures = c();
  temperaturesSRK = c();
  fractions = sapply(Substances, function(frac) frac$Fraction);
  Ac = sapply(Substances, function(ac) ac$Ac);
  Pc = sapply(Substances, function(pc) pc$Pc);
  Tc = sapply(Substances, function(tc) tc$Tc);
  substances = names(Substances)
  x = fractions;
  t0 = calc.bubble.T();
  Antoine.limit(t0, verbose=T);
  t0SRK = calc.bubbleSRK.T();
  SRK.limit(t0);
  dewpoint = calc.dewpoint.T();
  t100 = dewpoint$temperature;
  Antoine.limit(t100, verbose=T);
  dewpointSRK = calc.dewpointSRK.T();
  t100SRK = dewpointSRK$temperature;
  SRK.limit(t100SRK);
  printf('Bubble = %3.3f, Dew = %3.3f\t(Antoine)\n',t0,t100);
  printf('Bubble = %3.3f, Dew = %3.3f\t(SRK)\n',t0SRK,t100SRK);
  # 3. VLE calculation from T0 to T100 at P
  b = fractions # original concentration -> solution vector
  l = 1 # liquid part
  temperature = seq(t0,t100,1);
  temperatureSRK = seq(t0SRK,t100SRK,1);
  steps = length(temperature);
  stepsSRK = length(temperatureSRK);
  Boilingline = list(vapratio=0,
                     temperature=t0);
  BoilinglineSRK = list(vapratio=0,
                        temperature=t0SRK);
  xList = list(x=fractions);
  xListSRK = list(x=fractions);
  for (i in 2:(steps-1)){ # Antoine
    K = PsAntoine(temperature[i]) / pressure; # Equilibrium constants
    l = optimize(f=objF, lower=0, upper=1, tol=e)[[1]]; # one dimensional optimization
    xList = rbind(xList, objF(l, X=T));
    vapratio = 100 - l * 100;
    Boilingline = rbind(Boilingline, c(vapratio,
                                       temperature[i]));
  }
  Boilingline = rbind(Boilingline, c(100,
                                     t100));
  xList = rbind(xList,  dewpoint$x);
  for (i in 2:(stepsSRK-1)){ # SRK
    K = PsSRK(temperatureSRK[i]) / pressure; # Equilibrium constants
    l = optimize(f=objF, lower=0, upper=1, tol=e)[[1]]; # one dimensional optimization
    xListSRK = rbind(xListSRK, objF(l, X=T));
    vapratio = 100 - l * 100;
    BoilinglineSRK = rbind(BoilinglineSRK, c(vapratio,
                                             temperatureSRK[i]));
  }
  BoilinglineSRK = rbind(BoilinglineSRK, c(100,
                                           t100SRK));
  xListSRK = rbind(xListSRK,  dewpointSRK$x);
  mint0 = min(c(t0, t0SRK));
  maxt100 = max(c(t100, t100SRK));
  subtitle = c('p =', pressure, 'Pa');
  pdf(file='FACE2-csbl_n.pdf');
  plot(c(0,100), c(mint0,maxt100),
       main='FACE#2\nclosed system',
       xlab='Evaporated fraction [%]',
       ylab='T [K]',
       sub=paste(subtitle, collapse=' '),
       type='n')
  grid(NULL,NULL, lty = 6);
  legend('topleft', legend=c(paste('Dewpoint SRK:', round(t100SRK,1), 'K'),
                             paste('Dewpoint Antoine:', round(t100,1), 'K'),
                             paste('Bubblepoint SRK:', round(t0SRK,1), 'K'),
                             paste('Bubblepoint Antoine:', round(t0,1), 'K')),
  bty='n', cex=0.8, pch=c(2,2,1,1));
  points(c(0,100), c(t0,t100), pch=1) # Antoine
  lines(Boilingline,lwd=2);
  points(c(0,100), c(t0SRK,t100SRK), pch=2) # SRK
  lines(BoilinglineSRK, lwd=2, lty=2);
  dev.off();    
  return(list(Antoine=Boilingline,
              SRK=BoilinglineSRK));
}