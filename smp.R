smp = function(Substances=c(),
               Boilingline=c(),
               Pressure=NULL,
               UNIFAC=FALSE,
               monitor=FALSE,
               plot=FALSE){
  # Minimize the differantial between simulated and target boiling-, and
  # dew temperatures, by finding the optimal fractions of the entered
  # substances.
  #
  # Arguments:                                            
  #   Substances:     substances, c('ethanol', 'n-decane',...)
  #   Boiling Point:  target boiling temperature
  #   Dew Point:      target dew temperatur
  #   Pressure:       system pressure                                                        
  #
  # Returns:
  #   Vector with the molar fractions of the substances, (x1,x2,...,xn).
  #
  # Load and define functions:
  source('calc.bubblepoint.T.R')
  source('calc.dewpoint.T.R')
  source('Antoine.P.R')
  source('Antoine.T.R')
  #source('gamp.R')
  source('calc.boiling.T.R')
  source('calc.boilingline.R')
  source('sub.check.R')
  printf = function(...) cat(sprintf(...));
  Ps = function(Temperature, ...) {# list with saturated pressures
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature)
  }
  objF = function(x, ...){ # ideal assumptions
    x = abs(x);
    x = x / sum(x);
    Fractions = x;
    if(abs(1 - sum(x)) > 1e-7) warning("\nThe sum of the fractions differ 1.");
    T0 = calc.bubblepoint.T(Substances, Fractions, Pressure)
    T100 = calc.dewpoint.T(Substances, Fractions, Pressure)
    break
    b = Fractions; # original concentration -> solution vector
    l = 1 # liquid part
    lList = NULL # initiate list for liquid fraction
    cList = NULL # list for iteration control
    xList = NULL # list with liquid fractions
    if(length(Boilingline) == 2){
      returnVal = sum(abs(T0 - Boilingline[1]), abs(T100 - Boilingline[2])); 
    }
    else{
      Trange = Boilingline[,1];
      steps = length(Trange);
      objF = function(l, ...){ # function for optimization
        xf = l * (1 - K) + K;  # x factors for linear equation system
        A = matrix(0, nos, nos);  # initiate matrix
        for (i in 1:nos) {
          A[i, i] <- xf[i]; # fill the matrix
        }
        x = solve(A, b);  # solve the system
        abs(abs(sum(x) - 1) - abs(sum(K * x) - 1));
      }
      lList = l;
      for (i in 2:(steps - 1)){
        c = 0;
        K = Ps(Trange[i]) / Pressure; # Equilibrium constants
        l = optimize(f=objF, lower=0, upper=1, tol=e)[[1]]; # one dimensional optimization
        l1 = l;
        if(abs(0.5 - l) < e){ # starts iteration if the optimization fails
          l = 1;
          repeat{
            xf = l * (1 - K) + K;  # x factors for linear equation system
            A = matrix(0, nos, nos);  # initiate matrix
            for (j in 1:nos) {
              A[j, j] = xf[j]; # fill the matrix
            }
            x = solve(A, b);  # solve the system
            if((abs(sum(x) - 1) < 0.001) && (abs(sum(K * x) - 1) < e) || (c > 9999)){
              break;
            } else if (abs(sum(K * x) < 1)){
              l = l + e;
            }
            else{
              l = l - 0.001;
            }
            c = c + 1;
          }  
        } 
        xList = rbind(xList, x);    
        lList = rbind(lList, l);
        cList = rbind(cList, c);
      }
      l = 0;
      lList = rbind(lList, l);
      colnames(xList) = Substances;
      TemLiquid = cbind(Trange, lList);
      TemLiquid[1,1] = T0;
      TemLiquid[steps,1] = T100;
      TemVapor = TemLiquid;
      TemVapor[,2] = 1 - TemVapor[,2];
      colnames(TemLiquid) = c('T [K]','% liquid');
      returnVal = sum(abs(TemVapor[,2] - Boilingline[,2]),
                      abs(TemVapor[1,1] - Boilingline[1,1]),
                      abs(TemVapor[steps,1]- Boilingline[steps,1]));
    }
    if(monitor){
      points(TemVapor[,2],TemVapor[,1]);
      print(returnVal);
    }
    result = returnVal;
    return(result);
  }
  objFU = function(x, verbose=F, ...){ # UNIFAC liquid phase model
    x = abs(x);
    x = x / sum(x);
    Fractions = x;
    if(abs(1 - sum(x)) > 1e-7) warning("\nThe sum of the fractions differ 1.");
    Uvalues = UNIFAC.gen(Substances) # load or generate UNIFAC values
    unu = Uvalues[[1]]
    aij = Uvalues[[2]]
    T0 = calc.bubblepoint.T(Substances, Fractions, Pressure, UNIFAC=TRUE)$Temperature;
    T100 = calc.dewpoint.T(Substances, Fractions, Pressure, UNIFAC=TRUE)$Temperature;
    if(verbose) printf('Bubble = %3.3f, Dew = %3.3f\n',T0,T100);
    Trange = Boilingline[,1];
    steps = length(Trange);
    boilline = c()
    for(i in 2:(steps - 1)){
      Pd = calc.dewpoint.P(Substances,Fractions,Trange[i])$Pressure;
      Pb = calc.bubblepoint.P(Substances,Fractions,Trange[i])$pressure;
      # steam check
      if((Pd < Pressure) && (Pressure < Pb)){
        if(verbose) printf('\nsteam');
        v = (Pb - Pressure) / (Pb - Pd); # first guess
        boilline = c(boilline, Trange[i], v);
        x = Fractions;
      }
      else if((Pressure > Pb) && (verbose)) printf('\nliquid');
      if ((Pressure > Pd) && (verbose)) printf('\nvapor');
    }
    if(verbose)printf('\n');
    boilline = matrix(boilline, ncol=2, byrow=T);
    boilline = rbind(c(T0,0), boilline);
    boilline = rbind(boilline, c(T100,1));
    returnVal = sum(abs(boilline[,2] - Boilingline[,2]),
                    abs(boilline[1,1] - Boilingline[1,1]),
                    abs(boilline[steps,1]- Boilingline[steps,1]));
    if(monitor){
      points(boilline[,2], boilline[,1]);
      print(returnVal);
    }
    result = returnVal;
    return(result);
  } # end of objFU

  if (length(Boilingline) == 2){
    Tmin = Boilingline[1];
    Tmax = Boilingline[2];
  }
  Tmin = min(Boilingline[,1]);
  Tmax = max(Boilingline[,1]);
  s = Substances;
  P = Pressure;
  nos = length(s);
  e = 1e-7; # accuracy
  sublist = NULL;
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i]);
    sublist = rbind(sublist, tmp[1,]);
    sublist[i,7] = calc.boiling.T(Substances[i], P);
  }
  differ1 = max(sublist[,7]) - Tmax;
  printf('\nThe destination dewpoint is %3.1f K under the boilingpoint of the highboiler.',differ1)
  if(differ1 < 0){
    stop('\nThe Dewpoint has a difference of ',round(differ1,2),
         ' K to the boilingpoint of the highboiler.\nAbort blending.')
  }
  differ2 = Tmin - min(sublist[,7]);
  printf('\nThe destination bubblepoint is %3.1f K above the boilingpoint of the lowboiler.\n\nstart blending ...\n',differ2)
  if(differ2 < 0){
    stop('\nThe Boilingpoint has a difference of ',round(differ2,2),
         ' K to the boilingpoint of the lowboiler.\nAbort blending.')
  }
  m = differ2 - differ1;
  sublist = sublist[order(sublist$T.boil),]; # sort substances by boiling temperature
  for (i in 1:nos)Substances[i] = strsplit(as.character(sublist$Substance[i]),' ')[[1]][1]; # new order
  View(sublist);
  if(monitor){
    plot(Boilingline[,2],
         Boilingline[,1],
         pch=8,main='FACE#5 Boilingline',
         ylab='T [K]',xlab='evaporated fraction [%]');
    legend(x=0.5, y=600, c('experiment','simulation'),
           pch=c(8,1), bty='n');
  }
  fractions = rep(1,nos);
  fractions[which.max(sublist[,7])-1] = 3;
  fractions[which.max(sublist[,7])] = 6;
  fractions = fractions / sum(fractions);
  ptm = proc.time();
  if(UNIFAC){
    result = optim(par=fractions, fn=objFU);
  }
  else{
    result = optim(par=fractions, fn=objF);  
  }
  duration = proc.time() - ptm;
  fractions = abs(result$par);
  fractions = fractions / sum(fractions);
  boillines=NULL;
  if(plot){ # plot result
    dev.off();
    plotF = function(x, ...){
      x = abs(x);
      x = x / sum(x);
      Fractions = x;
      if(abs(1 - sum(x)) > 1e-7) warning("\nThe sum of the fractions differ 1.");
      T0 = calc.bubblepoint.T(Substances, Fractions, Pressure)
      T100 = calc.dewpoint.T(Substances, Fractions, Pressure)
      b = Fractions; # original concentration -> solution vector
      l = 1 # liquid part
      lList = NULL # initiate list for liquid fraction
      cList = NULL # list for iteration control
      xList = NULL # list with liquid fractions
      if(length(Boilingline) == 2){
        returnVal = sum(abs(T0 - Boilingline[1]), abs(T100 - Boilingline[2])); 
      }
      else{
        Trange = Boilingline[,1];
        steps = length(Trange);
        objF = function(l, ...){ # function for optimization
          xf = l * (1 - K) + K;  # x factors for linear equation system
          A = matrix(0, nos, nos);  # initiate matrix
          for (i in 1:nos) {
            A[i, i] <- xf[i]; # fill the matrix
          }
          x = solve(A, b);  # solve the system
          abs(abs(sum(x) - 1) - abs(sum(K * x) - 1));
        }
        lList = l;
        for (i in 2:(steps - 1)){
          c = 0;
          K = Ps(Trange[i]) / Pressure; # Equilibrium constants
          l = optimize(f=objF, lower=0, upper=1, tol=e)[[1]]; # one dimensional optimization
          l1 = l;
          if(abs(0.5 - l) < e){ # starts iteration if the optimization fails
            l = 1;
            repeat{
              xf = l * (1 - K) + K;  # x factors for linear equation system
              A = matrix(0, nos, nos);  # initiate matrix
              for (j in 1:nos) {
                A[j, j] = xf[j]; # fill the matrix
              }
              x = solve(A, b);  # solve the system
              if((abs(sum(x) - 1) < 0.001) && (abs(sum(K * x) - 1) < e) || (c > 9999)){
                break;
              } else if (abs(sum(K * x) < 1)){
                l = l + e;
              }
              else{
                l = l - 0.001;
              }
              c = c + 1;
            }  
          } 
          xList = rbind(xList, x);    
          lList = rbind(lList, l);
          cList = rbind(cList, c);
        }
        l = 0;
        lList = rbind(lList, l);
        colnames(xList) = Substances;
        TemLiquid = cbind(Trange, lList);
        TemLiquid[1,1] = T0;
        TemLiquid[steps,1] = T100;
        TemVapor = TemLiquid;
        TemVapor[,2] = 1 - TemVapor[,2];
        colnames(TemLiquid) = c('T [K]','% liquid');
      }
      result = list(boilExp=Boilingline,
                    boilSim=TemVapor);
      return(result);
    }
    boillines = plotF(fractions);
    pdf(file='Test1.pdf');
    plot(boillines$boilExp[,2],
         boillines$boilExp[,1],
         pch=8,main='FACE#5 Boilingline',
         ylab='T [K]',xlab='evaporated fraction [%]');
    points(boillines$boilSim[,2],
           boillines$boilSim[,1]);
    lines(boillines$boilSim[,2],
          boillines$boilSim[,1],
          lty=1, lwd=2);
    legend('topright', c('experiment','simulation'),
           lty=c(0,1), pch=c(8,1), bty='n');
    dev.off();
  }
  result = list(substances=Substances,
                fractions=fractions,
                boillines=boillines,
                precision=result$value,
                time=duration[3])
  # plot
#  calc.boilingline(Substances, fractions, Pressure, verbose=T);
  return(result);
}
#  result = optim(par=rep(1,nos)/nos, fn=objF, method='L-BFGS-B', lower=1e-2, upper=Inf);