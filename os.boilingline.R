os.boilingline = function(substances=NULL,
                          fractions=NULL,
                          pressure=NULL,
                          mixname=NULL,
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
  if (is.null(mixname)){
    warning('No name specified, named "unknown"');
    mixname='unknown'
  }
  subtitle = c('p =', pressure, 'Pa'); # subtitle for plotting 
  Substances = rep(list(''), nos); # initiate the list of lists
  for(i in 1:nos){ # collect substances properties
    Substances[[i]] = Substance(substances[i]); # get properties
    Substances[[i]]$MassFraction = fractions[i];
    Substances[[i]]$Mass = fractions[i] * 100;
    Substances[[i]]$Fraction = as.numeric(fractions[i] /
                                            Substances[[i]]$MolareMass); # temporary value
    A = Substances[[i]]$Antoine$A[1];
    B = Substances[[i]]$Antoine$B[1];
    C = Substances[[i]]$Antoine$C[1];
    tsat = Antoine.T(A, B, C, pressure); # calc boiling temperature using Antoine coefficients
    Substances[[i]]$Tsat = tsat;
    Ac = Substances[[i]]$Ac;
    Pc = Substances[[i]]$Pc;
    Tc = Substances[[i]]$Tc;
    TsatSRK = SRK(pressure, x=1, Ac=Ac, Pc=Pc, Tc=Tc, method = 'ts'); # calc boiling using cubic equation
    Substances[[i]]$TsatSRK = TsatSRK$temperature; 
  }
  abc=sapply(Substances, function(abc) abc$Antoine[1,1:3]);
  A = as.numeric(abc[1,]);
  B = as.numeric(abc[2,]);
  C = as.numeric(abc[3,]);
  names(Substances) = substances; # add the names
  Substances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  m = sum(sapply(Substances, function(mass) mass$Fraction)); 
  masses = sapply(Substances, function(mass) mass$Mass);
  mass = sum(masses)
  for(i in 1:nos){
    Substances[[i]]$MolareAmount = as.numeric(Substances[[i]]$Fraction / m);
    Substances[[i]]$Fraction = as.numeric(Substances[[i]]$MolareAmount);  
  }
  Suborg = Substances;
  Ac = sapply(Substances, function(ac) ac$Ac);
  Pc = sapply(Substances, function(pc) pc$Pc);
  Tc = sapply(Substances, function(tc) tc$Tc);
  substances = names(Substances)
  u = UNIFAC.gen(substances); # load UNIFAC values
  unu = u[[1]];
  aij = u[[2]];
  filename=paste(c(mixname,'-mass-progression.pdf'), collapse='');
  pdf(file=filename);
  par(mfrow=c(3,2), oma=c(4.5, 4, 4, 2.5), mar=rep(.1, 4), cex=0.7, las=1)
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
        result = SRK(pressure=pressure, x=x, Tc=Tc, Pc=Pc, Ac=Ac, method = 'bubbleT');
      } 
      if(r == 4){
        result = calc.bubbleSRK.T();
      }
      if(r == 5){
        result = calc.bubbleSRKU.T();
      }
      temperature = result$temperature;
      y = result$y; # molare vapor fractions
      Mm = as.numeric(sapply(Substances, function(molmass) molmass$MolareMass)); # get molare masses
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
      if(mass < 10){ # last points controll
        temp = length(temperatures);
        if(temperature < temperatures[temp]){ # check for raising temperature
          warning('Model#', r, ' temperature from ', round(temperature,2), 'K to ',
                  round(temperatures[temp],2), 'K corrected.');
          temperature = temperatures[temp]; # workaround to prevent something happen wich don't might happen
        }
        temp = length(vapor);
        if (round(vapor[temp], 0) == round(100 - mass, 0)){
          warning('Model#', r, ' vapor from ', round(100 - mass,0), '% to ',
                  round(101 - mass, 0),'% corrected.');
          mass = mass - 1;
        }
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
    if(r == 4) BlSRKpure = list(vapor=vapor,
                                temperature=temperatures);
    if(r == 5) BlSRKU = list(vapor=vapor,
                             temperature=temperatures);
    for(i in (nos - 1):1){
      massstack[i,] = massstack[(i + 1),] + massstack[i,]; # sum masses -> build the stack
    }
    grays = gray.colors(nos);
    if((r %% 2) && (r != 5)){ # alternating plotting
      plot(x = temperatures, y = massstack[1,],
           type='n',
           ann=FALSE,
           xaxt='n');
      legend('topright',
             legend=paste(c('model#',r),collapse = ''),
             bty='n',
             cex=0.9);
    }
    if(!(r %% 2)){
      plot(x = temperatures, y = massstack[1,],
           type='n',
           ann=FALSE,
           xaxt='n',
           yaxt='n');
      legend('topright',
             legend=paste(c('model#',r),collapse = ''),
             bty='n',
             cex=0.9);
    }
    if(r == 5){
      plot(x = temperatures, y = massstack[1,], # empty bottom left plot
           type='n',
           ann=FALSE);
      legend('topleft', 
             legend=substances,
             col=grays,
             lwd=2,
             bty='n', cex=0.8);
      legend('topright',
             legend=c('model#1 Antoine pure',
                      'model#2 Antoine UNIFAC',
                      'model#3 SRK',
                      'model#4 SRK pure',
                      'model#5 SRK UNIFAC'),
             bty='n',
             cex=0.8);
      axis(1, at=c(round(temperatures[1],0),round(Substances[[1]]$Tsat,0)));
      plot(x = temperatures, y = massstack[1,],
           type='n',
           ann=FALSE,
           yaxt='n');
      legend('topright',
             legend=paste(c('model#',r),collapse = ''),
             bty='n',
             cex=0.9);
      axis(1, at=c(round(temperatures[1],0),round(Substances[[1]]$Tsat,0)));
    }
    for(i in 1:nos){
      lines(x = temperatures, y = massstack[i,], col=grays[i]);
      xx = c(temperatures, rev(temperatures));
      yy = c(massstack[i,], rep(0,101));
      polygon(xx, yy, col=grays[i], border=NA);
    }
    Substances = Suborg; # reset substances list
    masses = sapply(Substances, function(mass) mass$Mass);
    mass = sum(masses);             
  } # end model loop -> start post processing
  title(paste(c(mixname,'\nProgression Curves'), collapse=''), outer=TRUE, cex=0.8)
  mtext(paste(c('T [K]\n',subtitle), collapse=' '), 1, 3, outer=TRUE,cex=0.8) # x-axis
  mtext('Substances progression [%]', 2, 3, outer=TRUE, las=0, cex=0.8) # y-axis
  dev.off();
  duration = proc.time() - ptm; # calculation time
  printf('\nCalculation takes %3.3f s\n', as.numeric(duration[1]));
  temp1 = length(BlAntoine$temperature);
  temp2 = length(BlUNIFAC$temperature);
  temp3 = length(BlSRK$temperature);
  temp4 = length(BlSRKpure$temperature);
  temp5 = length(BlSRKU$temperature);
  mint0 = min(c(BlAntoine$temperature[1],
                BlUNIFAC$temperature[1],
                BlSRK$temperature[1],
                BlSRKpure$temperature[1],
                BlSRKU$temperature[1]));
  maxt100 = max(c(BlAntoine$temperature[temp1],
                  BlUNIFAC$temperature[temp2],
                  BlSRK$temperature[temp3],
                  BlSRKpure$temperature[temp4],
                  BlSRKU$temperature[temp5]));
  pdf(file=paste(c(mixname,'-osbl.pdf'),collapse=''));
  plot(c(0,100), c(mint0, maxt100),
       main=paste(c(mixname,'\nDistillation Curves'),collapse = ''),
       xlab='Evaporated fraction [%]',
       ylab='T [K]',
       sub=paste(subtitle, collapse=' '),
       type='n');
  grid(NULL,NULL, lty = 6);
  legend('topleft', legend=c('Antoine pure',
                             'Antoine UNIFAC',
                             'SRK',
                             'SRK pure',
                             'SRK UNIFAC'),
         lty=c(1,2,6,4,8),
         col=c('black', 'red', 'darkblue', 'darkgreen', 'maroon'),
         lwd=2,
         bty='n', cex=0.9)
  lines(BlAntoine$vapor, BlAntoine$temperature, lwd=3);
  lines(BlUNIFAC$vapor, BlUNIFAC$temperature, lty=2, lwd=3, col='red');
  lines(BlSRK$vapor, BlSRK$temperature, lty=6, lwd=3, col='darkblue');
  lines(BlSRKpure$vapor, BlSRKpure$temperature, lty=8, lwd=3, col='darkgreen');
  lines(BlSRKU$vapor, BlSRKU$temperature, lty=4, lwd=3, col='maroon');
  dev.off();
  result = list(Antoine=BlAntoine,
                UNIFAC=BlUNIFAC,
                SRK=BlSRK,
                SRKpure=BlSRKpure,
                SRKUNIFAC=BlSRKU);
  return(result);
}