calc.boilingline = function(Substances=c(),
                            Fractions=c(),
                            Pressure=NULL,
                            UNIFAC=FALSE,
                            verbose=FALSE){
  # calculate and plot the boilingline of multicomponent mixtures
  # load and define functions:
  source('calc.dewpoint.T.R');
  source('calc.bubblepoint.T.R');
  source('sub.check.R');
  source('Antoine.P.R');
  Percent = function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  Ps = function(Temperature, ...) {# list with saturated pressures
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature)
  }
  printf = function(...) cat(sprintf(...));
  e <- 1e-7  # accuracy
  nos = length(Substances); # count number of substances
  if(is.null(Pressure)){
    stop('No pressure [Pa] specified. Abort.');
  }
  if (nos == 0){
    stop('No substance specified. Abort.');
  }
  if (length(Fractions) != nos){ # check: N substances == N fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(1 - sum(Fractions)) > e){
    stop('Sum of mole fractions not equal to 1. Abort.');
  }
  sublist = NULL
  for (i in 1:length(Substances)){
    tmp = sub.check(Substances[i])
    sublist = rbind(sublist, tmp[1,])
  }
  # non ideal solution
  if(UNIFAC){
    source('calc.dewpoint.P.R');
    source('calc.bubblepoint.P.R');
    source('UNIFAC.gen.R');
    source('UNIFAC.R');    
    Uvalues = UNIFAC.gen(Substances) # load or generate UNIFAC values
    unu = Uvalues[[1]]
    aij = Uvalues[[2]]
    T0 = calc.bubblepoint.T(Substances, Fractions, Pressure, UNIFAC=TRUE)[[1]]
    # 2. run simulation to calculate dewpoint temperature
    T100 = calc.dewpoint.T(Substances, Fractions, Pressure, UNIFAC=TRUE)[[1]]
    if(verbose) printf('Bubble = %3.3f, Dew = %3.3f\n',T0,T100);
    Trange = seq(T0,T100,1)
    steps = length(Trange)
    boilline = c()
    for(i in 1:steps){
      Pd = calc.dewpoint.P(Substances,Fractions,Trange[i])$Pressure;
      Pb = calc.bubblepoint.P(Substances,Fractions,Trange[i]);
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
    boilline[steps,] = c(T100,1);
    colnames(boilline) = c('T [K]', '% vapor');
    maintitle = paste(c(Substances, 'boilingline.unifac'), collapse='-');
    subtitle = c('P =', Pressure, 'Pa');
    filename = paste(c(maintitle,'pdf'),collapse='.');
    pdf(file=filename);
    plot(c(0,1), c(T0,T100),
         main=toupper(paste(maintitle, collapse='-')),
         xlab='Evaporated fraction',ylab='T [K]',
         sub=paste(subtitle, collapse=' '),
         type='p');
    lines(boilline[,2],boilline[,1],type='l',lwd=2);
    dev.off();
    T0100 = c(T0,T100); # T.bubble, T.dew
    result=list(BubbleAndDewpoint=T0100,
                TemperatureToVapor=boilline);
  }
  # ideal solution
  else{
    # 1. run simulation to calculate bubblepoint temperature
    T0 = calc.bubblepoint.T(Substances, Fractions, Pressure)
    # 2. run simulation to calculate dewpoint temperature
    T100 = calc.dewpoint.T(Substances, Fractions, Pressure)
    if(verbose) printf('Bubble = %3.3f, Dew = %3.3f\n',T0,T100);
    # 3. VLE calculation from T0 to T100 at P
    b = Fractions # original concentration -> solution vector
    l = 1 # liquid part
    lList = NULL # initiate list for liquid fraction
    cList = NULL # list for iteration control
    xList = NULL # list with liquid fractions
    Trange = seq(T0,T100,1)
    steps = length(Trange)
    objF = function(l, ...){ # function for optimization
      xf = l * (1 - K) + K;  # x factors for linear equation system
      A = matrix(0, nos, nos);  # initiate matrix
      for (i in 1:nos) {
        A[i, i] <- xf[i]; # fill the matrix
      }
      x = solve(A, b);  # solve the system
      abs(abs(sum(x) - 1) - abs(sum(K * x) - 1));
    }
    for (i in 1:steps){
      c = 0;
      K = Ps(Trange[i]) / Pressure; # Equilibrium constants
      l = optimize(f=objF, lower=0, upper=1, tol=e)[[1]]; # one dimensional optimization
      l1 = l;
      if(abs(0.5 - l) < 0.001){ # starts iteration if the optimization fails
        l = 1;
        repeat{
          xf = l * (1 - K) + K;  # x factors for linear equation system
          A = matrix(0, nos, nos);  # initiate matrix
          for (j in 1:nos) {
            A[j, j] = xf[j]; # fill the matrix
          }
          x = solve(A, b);  # solve the system
          if((abs(sum(x) - 1) < e) && (abs(sum(K * x) - 1) < e) || (c > 9999)){
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
      if(verbose){
        printf('\nl1 = %.3f, l2 = %.3f, iterations = %1.0f, progress = %s',
               l1, l, c, Percent(i / steps));        
      }   
      xList = rbind(xList, x);    
      lList = rbind(lList, l);
      cList = rbind(cList, c);
    }
    if(verbose)printf('\n');
    colnames(xList) = Substances;
    temLiquid = cbind(Trange, lList);
    temLiquid[steps, ] = c(T100,0);
#    temLiquid = rbind(temLiquid,c(T100,0));
    colnames(temLiquid) = c('T [K]','% liquid');
    T0100 = c(T0,T100); # T.bubble, T.dew
    # plot
    Substances='Test1'
    maintitle = paste(c(Substances, 'boilingcurve'), collapse='-');
    subtitle = c('P =', Pressure, 'Pa');
    filename = paste(c(maintitle,'pdf'),collapse='.');
    pdf(file=filename)
    plot(c(0,1), c(400,700),
#          main=toupper(paste(maintitle, collapse='-')),
           main='FACE#5',
         xlab='Evaporated fraction',ylab='T [K]',
         sub=paste(subtitle, collapse=' '),
         type='p')
    lines(1-temLiquid[,2],temLiquid[,1],type='l',lwd=2);
    dev.off();    
    result=list(BubbleAndDewpoint=T0100,
                LiquidFractions=xList,
                Iterations=cList,
                TemperatureToVapor=temLiquid);
  }
  return(result)
}
