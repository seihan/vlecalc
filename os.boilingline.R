os.boilingline = function(substances=NULL,
                          fractions=NULL,
                          pressure=NULL,
                          verbose=F){
  source('Substance.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('UNIFAC.R');
  source('UNIFAC.gen.R');
  source('SRK.R');
  Percent = function(x, digits = 0, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  PsAntoine = function(temperature, ...) {# list with saturated pressures
    ps = Antoine.P(A, B, C, temperature);
    return(ps);
  }
  TsAntoine = function(pressure, ...){
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
  calc.bubble.T = function(...){
    ts = TsAntoine(pressure);
    temperature = sum(x * ts);
    c = 0;
    while(c < 100){
      ps = PsAntoine(temperature);
      p = sum(x * ps);
      y = abs(x * ps / pressure); 
      rat = p / pressure; 
      if(abs(1 - rat) < e){break;}
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
  calc.bubbleUNIFAC.T = function(...){
    ts = TsAntoine(pressure);
    temperature = sum(x * ts);
    c = 0;
    while(c < 100){
      act = UNIFAC(Fractions=x, unu=unu, aij=aij, Temperature=temperature);
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
  calc.bubbleSRKU.T = function(...){
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
    if(c == 100) printf('Warning: Iteration limit (cubic). Calculated pressure = %3.3f Pa.\n', p);
    result = list(temperature=temperature,
                  y=abs(y));
    return(result);
  }
  calc.bubbleSRK.T = function(...){
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
    if(c == 100) printf('Warning: Iteration limit (cubic pure). Calculated pressure = %3.3f Pa.\n', p);
    result = list(temperature=temperature,
                  y=abs(y));
    return(result);
  }
  calc.bubble.hetero.T = function(...){
    ts = TsAntoine(pressure);
    temperature = sum(x * ts);
    phi = rep(1, nos);
    p = pressure;
    c = 1;
    while(c < 400){
      pformer = p;
      act = UNIFAC(Fractions=x, unu=unu, aij=aij, Temperature=temperature);
      ps = PsAntoine(temperature);
      p = sum(phi * act * x * ps);
      y = phi * act * x * ps / p;
      phi = SRK(pressure=p, temperature=temperature, x=x, Tc=Tc, Pc=Pc, Ac=Ac, hetero=T);
      rat = p / pressure; 
      diff = (abs(p - pformer));
      if ((c != 1) && (abs(1 - rat) < 1e-4)){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - rat) / rat;
        temperature = temperature + temp;
      }
    }
    if(c == 400) printf('Warning: Iteration limit (hetero). Calculated pressure = %3.3f Pa.\n', p)
    result = list(temperature=temperature,
                  y=y,
                  iterations=c);
    return(result);
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
  subtitle = c('p =', pressure, 'Pa'); # subtitle for plotting 
  Substances = rep(list(''), nos); # initiate the list of lists
  for(i in 1:nos){ # collect substances properties
    Substances[[i]] = Substance(substances[i]); # get properties
    Substances[[i]]$MassFraction = fractions[i];
    Substances[[i]]$Mass = fractions[i] * 100;
    Substances[[i]]$Fraction = as.numeric(fractions[i] /
                                            Substances[[i]]$MolarMass); # temporary value
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
  abc=sapply(Substances, function(abc) abc$Antoine[1,1:3]);
  A = as.numeric(abc[1,]);
  B = as.numeric(abc[2,]);
  C = as.numeric(abc[3,]);
  names(Substances) = substances;
  Substances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  m = sum(sapply(Substances, function(mass) mass$Fraction)); 
  masses = sapply(Substances, function(mass) mass$Mass);
  mass = sum(masses)
  for(i in 1:nos){
    Substances[[i]]$MolarAmount = as.numeric(Substances[[i]]$Fraction / m);
    Substances[[i]]$Fraction = as.numeric(Substances[[i]]$MolarAmount);  
  }
  Suborg = Substances;
  Ac = sapply(Substances, function(ac) ac$Ac);
  Pc = sapply(Substances, function(pc) pc$Pc);
  Tc = sapply(Substances, function(tc) tc$Tc);
  substances = names(Substances)
  u = UNIFAC.gen(substances); # load UNIFAC values
  unu = u[[1]];
  aij = u[[2]];
  for(r in 1:5){ # start model loop
    fractions = as.numeric(sapply(Substances, function(frac) frac$Fraction));
    x = fractions;
    temperatures = c(); # initiate temperature table
    vapor = 0; # initiate vapor table
    massstack = masses; # initiate masses table
    printf('\n\nModel %d of 5\n\nProgress: ',r);
    while(round(mass, 0) > 0){ # simulation start from 100 to 0g substances
      if(r == 1){
        result = calc.bubble.T();
      }
      if(r == 2){
        result = calc.bubbleUNIFAC.T();
      }
      if(r == 3){
        result = SRK(pressure=pressure, x=x, Tc=Tc, Pc=Pc, Ac=Ac);
      } 
      if(r == 4){
        result = calc.bubbleSRKU.T();
      }
      if(r == 5){
        result = calc.bubbleSRK.T();
      }
      temperature = result$temperature;
      if(mass < 10){
        temp = length(temperatures);
        if(temperature < temperatures[temp]){ # check for raising temperature
         temperature = temperatures[temp]; # workaround to prevent something happen wich don't might happen
         printf("\ncorrected!");
        }
      }
      y = result$y; # molare vapor fractions
      Mm = as.numeric(sapply(Substances, function(molmass) molmass$MolarMass)); # get molare masses
      ym = abs(y * Mm / sum(y * Mm)); # convert molare vapor to mass fraction
      for(i in 1:nos){
        Substances[[i]]$Mass = round(abs(as.numeric(Substances[[i]]$Mass - ym[i])), 9); # substract 1g vapor
      }
      masses = sapply(Substances, function(mass) mass$Mass); # get masses
      mass = sum(masses);
      xm = as.numeric(masses / mass); # calculate mass fractions
      x = xm / Mm / sum(xm / Mm); # convert mass to molare fractions
      for(i in 1:nos){
        Substances[[i]]$Fraction = x[i]; # store new molare fractions
      }
      temperatures = c(temperatures, temperature); # store temperature
      vapor = c(vapor, 100 - mass); # store vapor
      massstack = cbind(massstack,masses); # store masses
      printf('%s ', Percent((100 - mass) * 1e-2));
    } # simulation end
    if((r == 1) || (r == 2)){
      temperatures = c(temperatures, Substances[[1]]$Tsat); # set last temperature (boiling highboiler)
    }
    else{
      temperatures = c(temperatures, Substances[[1]]$TsatSRK);      
    }
    if(r == 1) BlAntoine = list(vapor=vapor,
                                temperature=temperatures);
    if(r == 2) BlUNIFAC = list(vapor=vapor,
                               temperature=temperatures);
    if(r == 3) BlSRK = list(vapor=vapor,
                            temperature=temperatures);
    if(r == 4) BlSRKU = list(vapor=vapor,
                             temperature=temperatures);
    if(r == 5) BlSRKpure = list(vapor=vapor,
                                temperature=temperatures);
    for(i in (nos - 1):1){
      massstack[i,] = massstack[(i + 1),] + massstack[i,]; # sum masses -> build the stack
    }
    filename=paste(c('FACE#1-osmasses-model-',r,'.pdf'), collapse='');
    grays = gray.colors(nos);
    pdf(file=filename);
    plot(x = temperatures, y = massstack[1,],
         main=paste(c('FACE#1\nopen system - model#', r), collapse=''),
         xlab='T [K]',
         ylab='Substances progression [%]',
         sub=paste(subtitle, collapse=' '),
         type='n');  
    grid(NULL,NULL, lty = 6);
    legend('topright', legend=substances,
           col=grays,
           lwd=2,
           bty='n', cex=0.9);
    for(i in 1:nos){
      lines(x = temperatures, y = massstack[i,], col=grays[i]);
      xx = c(temperatures, rev(temperatures));
      yy = c(massstack[i,], rep(0,101));
      polygon(xx, yy, col=grays[i], border=NA);
    }
    dev.off();
    Substances = Suborg; # reset substances list
    masses = sapply(Substances, function(mass) mass$Mass);
    mass = sum(masses);             
  } # end model loop -> start post processing
  duration = proc.time() - ptm; # calculation time
  printf('\nCalculation takes %3.3f s\n', as.numeric(duration[1]));
  temp1 = length(BlAntoine$temperature);
  temp2 = length(BlUNIFAC$temperature);
  temp3 = length(BlSRK$temperature);
  temp4 = length(BlSRKU$temperature);
  temp5 = length(BlSRKpure$temperature);
  mint0 = min(c(BlAntoine$temperature[1],
                BlUNIFAC$temperature[1],
                BlSRK$temperature[1],
                BlSRKU$temperature[1],
                BlSRKpure$temperature[1]));
  maxt100 = max(c(BlAntoine$temperature[temp1],
                  BlUNIFAC$temperature[temp2],
                  BlSRK$temperature[temp3],
                  BlSRKU$temperature[temp4],
                  BlSRKpure$temperature[temp5]));
  pdf(file='FACE1-osbl.pdf');
  plot(c(0,100), c(mint0, maxt100),
       main='FACE#1\nopen system',
       xlab='Evaporated fraction [%]',
       ylab='T [K]',
       sub=paste(subtitle, collapse=' '),
       type='n');
  grid(NULL,NULL, lty = 6);
  legend('topleft', legend=c('Antoine',
                             'Antoine UNIFAC',
                             'SRK',
                             'SRK UNIFAC',
                             'SRK pure'),
         lty=c(1,2,6,4,8),
         col=c('black', 'gray8', 'darkblue', 'darkgreen', 'maroon'),
         lwd=2,
         bty='n', cex=0.9)
  lines(BlAntoine$vapor, BlAntoine$temperature, lwd=3);
  lines(BlUNIFAC$vapor, BlUNIFAC$temperature, lty=2, lwd=3, col='red');
  lines(BlSRK$vapor, BlSRK$temperature, lty=6, lwd=3, col='darkblue');
  lines(BlSRKU$vapor, BlSRKU$temperature, lty=4, lwd=3, col='darkgreen');
  lines(BlSRKpure$vapor, BlSRKpure$temperature, lty=8, lwd=3, col='maroon');
  dev.off();
  result = list(Antoine=BlAntoine,
                UNIFAC=BlUNIFAC,
                SRK=BlSRK,
                SRKUNIFAC=BlSRKU,
                SRKpure=BlSRKpure);
  return(result);
}