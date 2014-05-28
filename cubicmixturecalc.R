# get the critical parameters for a substance from list
lookupcritical<-function(substance)
{
  substance <- substances[substance,]
  return(substance)
}

# estimate of Psat given T
Pestimate <- function (temperature,Tc,Pc,Wc)
{
  Temp <- log(Pc)
  Temp1 <- log(10)*(1-Tc/temperature)*(7+7*Wc)/3
  pressure <- exp(Temp+Temp1)
  return (pressure)
}

# estimate of Tsat given P
Testimate <- function (pressure,Tc,Pc,Wc)
{
  Temp <- log(pressure/Pc)
  Temp1 <- 1 - Temp*3/(log(10)*(7+7*Wc))
  temperature <- Tc/Temp1
  return(temperature)
}

# short-cut estimate of K values used for first iteration
Kestimate <- function (temperature,pressure,Tc,Pc,Wc)
{
  Temp <- log(Pc/pressure)
  Temp1 <- log(10)*(1-Tc/temperature)*(7+7*Wc)/3
  K <- exp(Temp+Temp1)
  return (K)
}
# calculate a and b parameters of pure components
calcab <- function(Tc,Pc,Wc,Tr,temperature,eos)
{
  if (eos==1)
  {
    kappa <- 0.480 + 1.574*Wc - 0.176*Wc*Wc
    alpha <- (1+kappa*(1-sqrt(Tr)))^2
    a <- 0.42748*R*R*Tc*Tc*alpha/(Pc*1e+5)
    b <- 0.08664*R*Tc/(Pc*1e+5)
    Xi <- -kappa*sqrt(Tr) / (1+kappa*(1-sqrt(Tr)))
  }
  if (eos==2)
  {
    kappa <- 0.37464 + 1.54226*Wc - 0.26992*Wc*Wc
    alpha <- (1+kappa*(1-sqrt(Tr)))^2
    a <- 0.45724*R*R*Tc*Tc*alpha/(Pc*1e+5)
    b <- 0.07780*R*Tc/(Pc*1e+5)
    Xi <- -kappa*sqrt(Tr) / (1+kappa*(1-sqrt(Tr)))
  }
  if (eos==3)
  {
    kappa <- 0.134 + 0.508*Wc - 0.0467*Wc*Wc
    alpha <- exp((2.00+0.836*Tr)*(1-(Tr)^kappa))
    a <- 0.45724*R*R*Tc*Tc*alpha/(Pc*1e+5)
    b <- 0.07780*R*Tc/(Pc*1e+5)
    Xi <- 0.836*Tr*(1-(kappa+1)*Tr^kappa - 2.00*kappa*Tr^kappa)
  }  
  return (c(a,b,Xi))
}

# calculate a and b parameters of mixture
calcabmix <- function (a11,b11,a22,b22,k12,x1)
{
   x2 <- 1 - x1
   a12 <- sqrt(a11*a22)*(1-k12)
   a <- x1*x1*a11 + 2*x1*x2*a12 +  x2*x2*a22
   b <- x1 * b11 + x2 * b22
  
  return (c(a,b,a12))
}
# calc a, b
absrk <- function (nkomp,temperature,x,Tc,Pc,Wc,kij)
{
  a <- 0
  b <- 0
  am <- c()
  alpha <- c()
  aij <- matrix(0,nkomp,nkomp)
  bi <- c()
#  i <- 1
  for (i in 1:nkomp)
  {
    print(i)
    am <- c(am, 0.48 + 1.57 * Wc[i] - 0.176 * (Wc[i])^2)
    alpha <- c(alpha, (1+am[i] * (1 - sqrt(temperature/Tc[i])))^2)
    aij[i,i] <- aij[i,i] + 0.42748*R^2*Tc[i]^2/(Pc[i]*1e+5*alpha[i])
    a <- a + x[i]^2 * aij[i,i]
    bi <- c(bi, 0.08664 * R * Tc[i]/(Pc[i]*1e+5))
    b <- b + x[i] * bi[i]
  }
print(aij)
  if (nkomp != 1)
  {
    i <- 1
    while (i < nkomp - 1)
    {
      j <- i + 1
      while (j < nkomp)
      {
        aij[i,i] <- sqrt(aij[i,i] * a[j,j] * (1 - kij[i,j]))
        aij[j,i] <- aij[i,j]
        a <- a + 2 * x[i] * x[j] * aij[i,j]
        j <- j + 1
      }
      i <- i + 1
    }
  }
  return (c(a,b,aij,bi,alpha,am))
}


# calculate coefficients in the cubic equation of state
calccoeffs <- function (A,B,eos)
{
  C2 <- 0
  C1 <- 0
  C0 <- 0
  if (eos==1)	
  {
    C2 <- -1
    C1 <- A - B - B*B
    C0 <- -A*B
  }
  else
  {
    C2 <- -1 + B
    C1 <- A - 3*B*B - 2*B
    C0 <- -A*B + B*B + B*B*B
  }
  
  return (c(C0,C1,C2))
}

# function to solve the cubic equation of state
cubic <- function (temperature,pressure,x1,critical1,critical2,k12,eos)
{
  pressure <- pressure * 1e+5
  Tc1 <- critical1[1]
  Pc1 <- critical1[2]
  Wc1 <- critical1[3]
  Tr1 <- critical1[4]
  Tc2 <- critical2[1]
  Pc2 <- critical2[2]
  Wc2 <- critical2[3]
  Tr2 <- critical2[4]
  phi1<-0
  phi2<-0
  phisecond1<-0
  phisecond2<-0
  
  # calculate the coefficients of the cubic equation 
  
  paramsab1 <- calcab(Tc1,Pc1,Wc1,Tr1,temperature,eos)
  paramsab2 <- calcab(Tc2,Pc2,Wc2,Tr2,temperature,eos)
  paramsab <- calcabmix(paramsab1[1],paramsab1[2],paramsab2[1],paramsab2[2],k12,x1)
  
  A <- paramsab[1]*pressure/(R*temperature)^2
  A11 <- paramsab1[1]*pressure/(R*temperature)^2
  A22 <- paramsab2[1]*pressure/(R*temperature)^2
  A12 <- paramsab[3]*pressure/(R*temperature)^2
  B <- paramsab[2]*pressure/(R*temperature)
  B11 <- paramsab1[2]*pressure/(R*temperature)
  B22 <- paramsab2[2]*pressure/(R*temperature)
  Xi1 <- paramsab1[3]
  Xi2 <- paramsab2[3]
  coeffs <- calccoeffs(A,B,eos)
  
  # solve the cubic equation for Z
  
  Z <- solvecubic(coeffs[1],coeffs[2],coeffs[3])
  
  # calculate fugacity coefficients
  
  if (Z[3]>0)
  {
    phi1 <- calcfugacity(Z[1],A,B,x1,B11,A11,A12,eos)
    phi2 <- calcfugacity(Z[1],A,B,x1,B22,A12,A22,eos)
    phisecond1 <- calcfugacity(Z[3],A,B,x1,B11,A11,A12,eos)
    phisecond2 <- calcfugacity(Z[3],A,B,x1,B22,A12,A22,eos)
  }
  else
  {
    phi1 <- calcfugacity(Z[1],A,B,x1,B11,A11,A12,eos)
    phi2 <- calcfugacity(Z[1],A,B,x1,B22,A12,A22,eos)
    Z[3]<-0
    phisecond1<-0 
    phisecond2<-0
  }
  
  # calculate enthalpy and entropy departure functions
  
  if (Z[1]>0)
  {
    Hdep0 <- calcenthalpy(temperature,Z[1],x1,A,B,A11,A22,A12,Xi1,Xi2,eos)
    Sdep0 <- calcentropy(temperature,Z[1],x1,A,B,A11,A22,A12,Xi1,Xi2,eos)
  }
  else
  {
    Hdep0 <- 0
    Sdep0 <- 0
  }
  
  if (Z[3]>0)
  {
    Hdep2 <- calcenthalpy(temperature,Z[3],x1,A,B,A11,A22,A12,Xi1,Xi2,eos)
    Sdep2 <- calcentropy(temperature,Z[3],x1,A,B,A11,A22,A12,Xi1,Xi2,eos)
  }
  else
  {
    Hdep2 <- 0
    Sdep2 <- 0
  }
  
  # return values to main routine
  
  return (c(Z[1],phi1,phi2,Z[3],phisecond1,phisecond2,Hdep0,Hdep2,Sdep0,Sdep2))
}
# function for solving cubic equations
solvecubic <- function (C0,C1,C2)
{
  Q1 <- C2*C1/6 - C0/2 - C2*C2*C2/27
  P1 <- C2*C2/9 - C1/3
  D <- Q1*Q1 - P1*P1*P1
  if (D>=0)
  {
    temp1 <- abs(Q1+sqrt(D))^0.3333333333
    temp1 <- temp1 * (Q1 + sqrt(D)) / abs(Q1 + sqrt(D))
    temp2 <- abs(Q1-sqrt(D))^0.3333333333
    temp2 <- temp2 * (Q1 - sqrt(D))/abs(Q1 - sqrt(D))
    Z0 <- temp1 + temp2 - C2/3
    Z1 <- 0
    Z2 <- 0
  }
  else
  {
    temp1 <- Q1*Q1/(P1*P1*P1)
    temp2 <- sqrt(1-temp1)/sqrt(temp1)
    temp2 <- temp2 * Q1/abs(Q1)
    Phi <- atan(temp2)
    if (Phi<0)
    {
      Phi<-Phi+pi
    }
    Z0 <- 2*sqrt(P1)*cos(Phi/3) - C2/3
    Z1 <- 2*sqrt(P1)*cos((Phi+2*pi)/3) - C2/3
    Z2 <- 2*sqrt(P1)*cos((Phi+4*pi)/3) - C2/3
    if (Z0<Z1)
    {
      temp <- Z0
      Z0 <- Z1
      Z1 <- temp
    }
    if (Z1<Z2)
    {
      temp <- Z1
      Z1 <- Z2
      Z2 <- temp
    }
    if (Z0<Z1)
    {
      temp <- Z0
      Z0 <- Z1
      Z1 <- temp
    }
  }
  return (c(Z0,Z1,Z2,D))
}

# calculate fugacity coefficient of species
calcfugacity <- function (Z,A,B,x1,Bi,A1i,A2i,eos)
{
  if (eos==1)
  {
     temp1 <- -log(Z-B) + Bi*(Z-1)/B
     temp3 <- log(1+B/Z)
     x2 <- 1-x1
     temp2 <- A*(2*(x1*A1i + x2*A2i)/A - Bi/B)/B
     phi <- exp(temp1-temp2*temp3)
  }
  else
  {
     temp1 <- -log(Z-B) + Bi*(Z-1)/B
     temp3 <- log((Z+2.41421356*B)/(Z-0.41421356*B))
     x2 <- 1-x1
     temp2 <- A*(2*(x1*A1i + x2*A2i)/A - Bi/B)/(2.82842713*B)
     phi <- exp(temp1-temp2*temp3)
  }
  
  return (phi)
}

# calculate enthalpy departure function
calcenthalpy <- function (temperature,Z,x1,A,B,A11,A22,A12,Xi1,Xi2,eos)
{
  x2 <- 1-x1
  Hdep <- 0
  temp3 <- x1*x1*A11*Xi1 + x2*x2*A22*Xi2 + x1*x2*A12*(Xi1+Xi2)
  
  if (eos==1)
  {
    temp4 <- log(1+B/Z)
    Hdep <- R*temperature*((Z-1) + (temp3-A)*temp4/B)
  }
  else
  {
    temp4 <- log((Z+2.41421356*B)/(Z-0.41421356*B))
    Hdep <- R*temperature*((Z-1) + (temp3-A)*temp4/(2.82842713*B))
  }
  
  return (Hdep)
}

# calculate entropy departure function
calcentropy <- function (temperature,Z,x1,A,B,A11,A22,A12,Xi1,Xi2,eos)
{
  x2 <- 1-x1
  Sdep <- 0
  temp3 <- x1*x1*A11*Xi1 + x2*x2*A22*Xi2 + x1*x2*A12*(Xi1+Xi2)
  
  if (eos==1)
  {
     temp4 <- log(1+B/Z)
    Sdep <- R*(log(Z-B) + temp3*temp4/B)
  }
  else
  {
     temp4 <- log((Z+2.41421356*B)/(Z-0.41421356*B))
    Sdep <- R*(log(Z-B) + temp3*temp4/(2.82842713*B))
  }
  
  return (Sdep)
}
# main routine to do the calculation
mixturecalc <- function (z1,substance1,substance2,k12,temperature,pressure,eos,calc)
{
  substance <- lookupcritical(substance1)
  Tc1 <- as.numeric(substance[2])
  Pc1 <- as.numeric(substance[3])
  Wc1 <- as.numeric(substance[4])
  Tr1 <- temperature/Tc1
  critical1 <- c(Tc1,Pc1,Wc1,Tr1)
  substance <- lookupcritical(substance2)
  Tc2 <- as.numeric(substance[2])
  Pc2 <- as.numeric(substance[3])
  Wc2 <- as.numeric(substance[4])
  Tr2 <- temperature/Tc2
  critical2 <- c(Tc2,Pc2,Wc2,Tr2)
  temp2 <- 0

  # [1] perform bubble point temperature calculation
  if (calc==1)
  {
    P <- pressure
    x1 <- z1
    x2 <- 1-x1
#     temperature <- readline('Enter initial guess for T (K). [Default is to use last value. Program will make an initial guess if you enter zero or leave blank] ')
#     temperature <- as.numeric(temperature)
#     if (is.na(temperature)||(temperature<=0))
#       {
        temp1 <- Testimate(P,Tc1,Pc1,Wc1)
        temp2 <- Testimate(P,Tc2,Pc2,Wc2)
        temperature <- x1*temp1 + x2*temp2
#         print(temperature)
    #  }
    firstguess <- temperature
    K1 <- Kestimate(temperature,P,Tc1,Pc1,Wc1)
    K2 <- Kestimate(temperature,P,Tc2,Pc2,Wc2)
    y1 <- K1*x1
    y2 <- K2*x2
    Tformer <- 0
    i <- 0
    repeat
    {
      Tformer <- temperature
      y1former <- y1
      i <- i+1
      results <- cubic(temperature,P,x1,critical1,critical2,k12,eos)
      if (results[4]>0)
      {
        Z2 <- results[4]
        phisecond1 <- results[5]
        phisecond2 <- results[6]
        Hdep2 <- results[8]
        Sdep2 <- results[10]
      }
      else
      {
        Z2 <- results[1]
        phisecond1 <- results[2]
        phisecond2 <- results[3]
        Hdep2 <- results[7]
        Sdep2 <- results[9]
      }
      
      for (j in 1:25)
      {
        temp1 <- y1+y2
        y1 <- y1/temp1
        y2 <- y2/temp2
        results <- cubic(temperature,P,y1,critical1,critical2,k12,eos)
        Z0 <- results[1]
        phi1 <- results[2]
        phi2 <- results[3]
        Hdep0 <- results[7]
        Sdep0 <- results[9]
        y1 <- (x1*phisecond1)/phi1
        y2 <- (x2*phisecond2)/phi2
      }
      temp1 <- y1+y2
      temp2 <- 0.1*temperature*(1-temp1)/temp1
      temperature <- temperature + temp2
      print(temperature)
      if(((abs(temperature-Tformer)<0.0001)||(abs(y1-y1former)<0.000001))){break}
      ##i>19)
      
    } 
    temp3 <- abs(y1+y2-1)
    temp1 <- abs(x1*phisecond1 - y1*phi1)
    temp2 <- abs(x2*phisecond2 - y2*phi2)
    if ((temp1>0.00001)||(temp2>0.00001))
    {
      print('Warning: VLE part not fully converged')
    }
    if (temp3>0.00005)
    {
      print('Warning: not fully converged with this initial T estimate')
    }
    if ((!is.numeric(temperature))||(abs(Z0-Z2)<0.00005))
    {
      print('Sorry: calculation did not converge using the initial estimate')
      temperature <- firstguess
    }
    liqfn <- 1
    xfeed1 <- x1
    return(c(temperature,Z0,Z2,phi1,phi2,phisecond1,phisecond2,i))
  }

  # [2] perform bubble point pressure calculation
  if (calc==2)
  {
    x1 <- z1
    x2 <- 1-x1
    P <- pressure
    temp1 <- Pestimate(temperature,Tc1,Pc1,Wc1)
    temp2 <- Pestimate(temperature,Tc2,Pc2,Wc2)
    P <- x1*temp1 + x2*temp2
    firstguess <- P
    K1 <- Kestimate(temperature,P,Tc1,Pc1,Wc1)
    K2 <- Kestimate(temperature,P,Tc2,Pc2,Wc2)
    y1 <- K1*x1
    y2 <- K2*x2
    Pformer <- 0
    i <- 0
    repeat
    {
      Pformer <- P
      y1former <- y1
      i <- i+1
      results <- cubic(temperature,P,x1,critical1,critical2,k12,eos)
      if (results[4]>0)
      {
        Z2 <- results[4]
        phisecond1 <- results[5]
        phisecond2 <- results[6]
        Hdep2 <- results[8]
        Sdep2 <- results[10]
      }
      else
      {
        Z2 <- results[1]
        phisecond1 <- results[2]
        phisecond2 <- results[3]
        Hdep2 <- results[7]
        Sdep2 <- results[9]
      }
      for (j in 1:25)
      {
        temp1 <- y1+y2
        y1 <- y1/temp1
        y2 <- y2/temp1
        results <- cubic(temperature,P,y1,critical1,critical2,k12,eos)
        Z0 <- results[1]
        phi1 <- results[2]
        phi2 <- results[3]
        Hdep0 <- results[7]
        Sdep0 <- results[9]
        y1 <- x1*phisecond1/phi1
        y2 <- x2*phisecond2/phi2
      }
      temp2 <- P*(y1+y2-1)
      P <- P+temp2
      if((abs(P/Pformer-1)<0.000005)||(abs(y1-y1former)<0.000001)||(i>19)){break}
      print(P)
    }
    temp3 <- abs(y1+y2-1)
    temp1 <- abs(x1*phisecond1 - y1*phi1)
    temp2 <- abs(x2*phisecond2 - y2*phi2)
    if ((temp1>0.00001)||(temp2>0.00001))
    {
      print('Warning: VLE part not fully converged')
    }
    if (temp3>0.00005)
    {
      print('Warning: not fully converged with this initial P estimate')
    }
    if ((!is.numeric(P))||(abs(Z0-Z2)<0.00005))
    {
      print('Sorry: calculation did not converge using the initial estimate')
      P <- firstguess
    }
    liqfn <- 1
    xfeed1 <- x1	
    return(P)
  }

  # [3] perform dew point temperature calculation
  if(calc == 3)
  {
    P <- pressure
    y1 <- z1
    y2 <- 1-y1
    temp1 <- Testimate(P,Tc1,Pc1,Wc1)
    temp2 <- Testimate(P,Tc2,Pc2,Wc2)
    temperature <- y1*temp1 + (1-y1)*temp2
    firstguess <- temperature
    K1 <- Kestimate(temperature,P,Tc1,Pc1,Wc1)
    K2 <- Kestimate(temperature,P,Tc2,Pc2,Wc2)
    x1 <- y1/K1
    x2 <- y2/K2
    Tformer <- 0
    i <- 0
    repeat
    {
      Tformer <- temperature
      x1former <- x1
      i <- i+1
      results <- cubic(temperature,P,y1,critical1,critical2,k12,eos)
      Z0 <- results[1]
      phi1 <- results[2]
      phi2 <- results[3]
      Hdep0 <- results[7]
      Sdep0 <- results[9]
      for (j in 1:24)
      {
        temp1 <- x1+x2
        x1 <- x1/temp1
        x2 <- x2/temp1
        results <- cubic(temperature,P,x1,critical1,critical2,k12,eos)
        if (results[4]>0)
        {
          Z2 <- results[4]
          phisecond1 <- results[5]
          phisecond2 <- results[6]
          Hdep2 <- results[8]
          Sdep2 <- results[10]
        }
        else
        {
          Z2 <- results[1]
          phisecond1 <- results[2]
          phisecond2 <- results[3]
          Hdep2 <- results[7]
          Sdep2 <- results[9]
        }
        x1 <- y1*phi1/phisecond1
        x2 <- y2*phi2/phisecond2
      }
      temp1 <- x1+x2
      temp2 <- 0.1*temperature*(temp1-1)/temp1
      temperature <- temperature + temp2
      if((abs(temperature-Tformer)<0.0001)||(abs(x1-x1former)<0.000001)||(i>19)){break}
    } 
    temp3 <- abs(x1+x2-1)
    temp1 <- abs(x1*phisecond1 - y1*phi1)
    temp2 <- abs(x2*phisecond2 - y2*phi2)
    if ((temp1>0.00001)||(temp2>0.00001))
    {
      print('Warning: VLE part not fully converged')
    }
    if (temp3>0.00005)
    {
      print('Warning: not fully converged with this initial T estimate')
    }
    if ((!is.numeric(temperature))||(abs(Z0-Z2)<0.00005))
    {
      print('Sorry: calculation did not converge using the initial estimate')
      temperature <- firstguess
    }
    liqfn <- 0
    xfeed1 <- y1
    return(c(temperature,i))
  }
 
  # [4] perform dew point pressure calculation
  if(calc == 4)
  {
    y1 <- z1
    y2 <- 1-y1
    P <- pressure
    temp1 <- Pestimate(temperature,Tc1,Pc1,Wc1)
    temp2 <- Pestimate(temperature,Tc2,Pc2,Wc2)
    P <- 1/(y1/temp1 + y2/temp2)
    firstguess <- P
    K1 <- Kestimate(temperature,P,Tc1,Pc1,Wc1)
    K2 <- Kestimate(temperature,P,Tc2,Pc2,Wc2)
    x1 <- y1/K1
    x2 <- y2/K2
    Pformer <- 0
    i <- 0
    repeat
    {
      Pformer <- P
      x1former <- x1
      i <- i+1
      results <- cubic(temperature,P,y1,critical1,critical2,k12,eos)
      Z0 <- results[1]
      phi1 <- results[2]
      phi2 <- results[3]
      Hdep0 <- results[7]
      Sdep0 <- results[9]
      for (j in 0:24)
      {
        temp1 <- x1+x2
        x1 <- x1/temp1
        x2 <- x2/temp2
        results <- cubic(temperature,P,x1,critical1,critical2,k12,eos)
        if (results[4]>0)
        {
          Z2 <- results[4]
          phisecond1 <- results[5]
          phisecond2 <- results[6]
          Hdep2 <- results[8]
          Sdep2 <- results[10]
        }
        else
        {
          Z2 <- results[1]
          phisecond1 <- results[2]
          phisecond2 <- results[3]
          Hdep2 <- results[7]
          Sdep2 <- results[9]
        }
        x1 <- y1*phi1/phisecond1
        x2 <- y2*phi2/phisecond2
      }
      temp2 <- P*(1-x1-x2)
      P <- P+temp2
      if((abs(P/Pformer-1)<0.000005)||(abs(x1-x1former)<0.000001)||(i>19)){break}
    }
    temp3 <- abs(x1+x2-1)
    temp1 <- abs(x1*phisecond1 - y1*phi1)
    temp2 <- abs(x2*phisecond2 - y2*phi2)
    if ((temp1>0.00001)||(temp2>0.00001))
    {
      print('Warning: VLE part not fully converged')
    }
    if (temp3>0.00005)
    {
      print('Warning: not fully converged with this initial P estimate')
    }
    if ((!is.numeric(P))||(abs(Z0-Z2)<0.00005))
    {
      print('Sorry: calculation did not converge using the initial estimate')
      P <- firstguess
    }
    liqfn <- 0
    xfeed1 <- y1
    return(c(P,i))
  }

  # [5] perform isothermal VLE flash calculation  
  if (calc==5)
  {
    #pressure <- round(Pestimate(temperature,Tc,Pc,Wc)/10)*10
    xfeed1 <- z1
    results <- cubic(temperature,pressure,xfeed1,critical1,critical2,k12,eos)
    Z0 <- results[1]
    phi1 <- results[2]
    phi2 <- results[3]
    Hdep0 <- results[7]
    Sdep0 <- results[9]
    Z2  <- results[4]
    phisecond1  <- results[5]
    phisecond2  <- results[6]
    Hdep2  <- results[8]
    Sdep2  <- results[10]
    if ((Z0!=0)&&(Z2!=0))
    {
      K1  <- Kestimate(temperature,pressure,Tc1,Pc1,Wc1)
      K2  <- Kestimate(temperature,pressure,Tc2,Pc2,Wc2)
      for (i in 1:200)
      {
        liqfn <- (xfeed1*(K1-K2)-K1+K1*K2)/((1-K1)*(1-K2))
        x1 <- xfeed1/(liqfn+K1*(1-liqfn))
        y1 <- K1*x1
        results <- cubic(temperature,pressure,y1,critical1,critical2,k12,eos)
        Z0 <- results[1]
        phi1 <- results[2]
        phi2 <- results[3]
        Hdep0 <- results[7]
        Sdep0 <- results[9]
        results <- cubic(temperature,pressure,x1,critical1,critical2,k12,eos)
        Z2 <- results[3]
        phisecond1 <- results[5]
        phisecond2 <- results[6]
        Hdep2 <- results[8]
        Sdep2 <- results[10]
        K1 <- K1*x1*phisecond1/(y1*phi1)
        K2 <- K2*(1-x1)*phisecond2/((1-y1)*phi2)
      }
      y2 <- 1-y1
      x2 <- 1-x1
      temp1 <- abs(x1*phisecond1-y1*phi1)
      temp2 <- abs(x2*phisecond2-y2*phi2)
      if ((temp1>0.00001)||(temp2>0.00001))
      {
        print('Sorry- calculation did not converge')
      }
      if (!is.numeric(liqfn))
      {
        print('Sorry- no flash calculation possible under these conditions')
      }	
      if ((liqfn<0)||(liqfn>1))
      {
        print('Sorry - flash calculation does not give valid liquid fraction under these conditions')
      }
    }
    else
    {
      print("feed is a single phase")
      if (Z0!=0)
      {
        y1 <- xfeed1
        y2 <- 1-y1
        x1 <- 0
        x2 <- 0
        liqfn <- 0
      }
      else
      {
        x1 <- xfeed1
        x2 <- 1 - x1
        y1 <- 0
        y2 <- 0
        liqfn <- 1
      }
    }
    return(c(x1,y1))
  }

#   # [6] perform isothermal VLE calculation  
#   if (calc==6)
#   {
#     x1 <- z1
#     x2 <- 1-x1
#     temp1 <-Pestimate(temperature,Tc1,Pc1,Wc1)
#     temp2 <-Pestimate(temperature,Tc2,Pc2,Wc2)
#     P <- x1*temp1 + x2*temp2
#     results <- cubic(temperature,P,x1,critical1,critical2,k12,eos)
#     Z0 <- results[1]
#     phi1 <- results[2]
#     phi2 <- results[3]
#     Hdep0 <- results[7]
#     Sdep0 <- results[9]
#     Z2  <- results[4]
#     phisecond1  <- results[5]
#     phisecond2  <- results[6]
#     Hdep2  <- results[8]
#     Sdep2  <- results[10]
#     K1  <- Kestimate(temperature,P,Tc1,Pc1,Wc1)
#     K2  <- Kestimate(temperature,P,Tc2,Pc2,Wc2)
#     S <- 2
#     e <- 0.00001
#     for (var i = 0; i<200; i++)
#     {
#       liqfn <- (xfeed1*(K1-K2)-K1+K1*K2)/((1-K1)*(1-K2))
#       x1 <- xfeed1/(liqfn+K1*(1-liqfn))
#       y1 <- K1*x1
#       results <- cubic(temperature,P,y1,critical1,critical2,k12,eos)
#       Z0 <- results[1]
#       phi1 <- results[2]
#       phi2 <- results[3]
#       Hdep0 <- results[7]
#       Sdep0 <- results[9]
#       results <- cubic(temperature,P,x1,critical1,critical2,k12,eos)
#       Z2 <- results[3]
#       phisecond1 <- results[5]
#       phisecond2 <- results[6]
#       Hdep2 <- results[8]
#       Sdep2 <- results[10]
#       K1 <- K1*x1*phisecond1/(y1*phi1)
#       K2 <- K2*(1-x1)*phisecond2/((1-y1)*phi2)
#     }
#    return(c(x1,y1,P,phi1,phi2))
#  }
}

# global environment
substances<-read.csv("~/Dropbox/MA-Thesis/Fuels/calculations/values/criticalvalues.csv",head=T,sep=",")
R <- 8.31433
#Trange <- seq(260,360,1)
#txy<-c()
#for (i in 1:length(Trange))
#{
#  txy<-c(txy, Trange[1], mixturecalc(0.8,19,20,0,Trange[1],1,1,5))
#}