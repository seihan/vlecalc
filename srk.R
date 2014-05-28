##########################################################
##       Soave Redlich Kwong Equation of State          ##
##                                                      ##
##   -----------------------------------------------    ##
##   Program calculates the vapor mole fraction from    ##
##   from given liquid fraction and temperature for     ##
##   binary mixtures.                                   ##
##   The cubic eos by Soave, Redlich and Kwong based    ##
##   on the genarally form:                             ##
##     z = R*T / v-b - a(T) / (v² + b*v) * v / R*T      ##
##   is used.                                           ##
##   Besides the functions for calculating the          ##
##   coefficents is "mixturecalc" the function of       ##
##   interest. It returns the gas-fraction and P        ##
##   by default. Inside are the fugacity coefficents    ##
##   and molare volume for the liquid- and gas-phase    ##
##   calculate, which can easily returned too.          ##
##   -----------------------------------------------    ##
##                                                      ##
##                 MIXTURECALC                          ##
##   input:  x,substance1,substance2,binary parameter,  ##
##           T                                          ##
##   output: y,P,(phi_liquid[i],phi_gasous[i],v_l,v_g)  ##
##                                                      ##
##########################################################

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

# calculate a and b parameters of pure components
calcab <- function(Tc,Pc,Wc,Tr,temperature)
{
  R <- R * 1e-2
  kappa <- 0.48 + 1.574*Wc - 0.176*Wc*Wc
  alpha <- (1+kappa*(1-sqrt(Tr)))^2
  a <- 0.42748*R*R*Tc*Tc*alpha/(Pc)
  b <- 0.08664*R*Tc/(Pc)
  Xi <- -kappa*sqrt(Tr) / (1+kappa*(1-sqrt(Tr)))
  return (c(a,b,Xi))
}

# calculate a and b parameters of mixture
calcabmix <- function (a11,b1,a22,b2,k12,x1)
{
  x2 <- 1 - x1
  a12 <- sqrt(a11*a22)*(1-k12)
  a <- x1*x1*a11 + 2*x1*x2*a12 +  x2*x2*a22
  b <- x1 * b1 + x2 * b2
  return (c(a,b,a12))
}
# calculate molare volume by solving the SRK eos [1]
volsrk <- function (temperature,P,a,b,phase)
{
  R <- R * 1e-2
  pstr <- -(P*b^2+R*temperature*b-a) / 3 / P - (R*temperature / 3 / P)^2
  qstr <- -(R*temperature/3/P)^3 - R*temperature * (P*b^2 + R*temperature*b - a) / 6 / P^2 - a*b / 2 / P
  diskr <- pstr^3 + qstr^2
  if (diskr < 0)
  {
    if (qstr < 0){rstr <- - sqrt(abs(pstr))}
    else{rstr <- sqrt(abs(pstr))}
    cosphi <- qstr / rstr^3
    phi <- acos(cosphi)
    x1 <- -2 * rstr * cos(phi/3)
    x2 <- 2 * rstr * cos((pi - phi) / 3)
    x3 <- 2 * rstr * cos((pi + phi) / 3)
    if (phase == -1){v <- max(x1,x2,x3)} # [g]
    if (phase == 1){v <- min(x1,x2,x3)} # [l]
  }
  else
  {
    h1 <- - qstr + sqrt(diskr)
    h2 <- - qstr - sqrt(diskr)
    if (h1 < 0){h3 <- -1}
    else{h3 <- 1}
    if (h2 < 0){h4 <- -1}
    else{h4 <- 1}
    v <- h3 * (abs(h1))^(1/3) + h4 * (abs(h2))^(1/3)
  }
  vol <- v + R*temperature / 3 / P
  return(vol)
}

# calculate the fugacity from gas- or liquid phase
phisrk <- function (temperature,P,v,a,b,a1i,a2i,bi,x1)
{
  R <- R * 1e-2
  x2 <- 1 - x1
  h1 <- log(v/(v-b))
  h2 <- (a*bi)/(R*temperature*b^2)*(log((v+b)/v)-b/(v+b))
  h3 <- (bi/(v-b))-log(P*v/(R*temperature))
  h4 <- -2*(x1*a1i+x2*a2i)/(R*temperature*b)*log((v+b)/v)
  phi<- h1+h2+h3+h4
  return (exp(phi))
}
# calculate y, P, phi*, v* for a binary mixture
mixturecalc <- function(x,substance1,substance2,k12,temperature)
{
  substance <- lookupcritical(substance1)
  Tc1 <- as.numeric(substance[2])
  Pc1 <-  as.numeric(substance[3])
  Wc1 <-  as.numeric(substance[4])
  Tr1 <- temperature / Tc1
  substance <- lookupcritical(substance2)
  Tc2 <- as.numeric(substance[2])
  Pc2 <-  as.numeric(substance[3])
  Wc2 <-  as.numeric(substance[4])
  Tr2 <- temperature / Tc2
  ab1 <- calcab(Tc1,Pc1,Wc1,Tr1,temperature)
  ab2 <- calcab(Tc2,Pc2,Wc2,Tr2,temperature)
  # Liquid
  abmix <- calcabmix(ab1[1],ab1[2],ab2[1],ab2[2],k12,x)
  al <- abmix[1]
  a11 <- ab1[1]
  al12 <- abmix[3]
  a22 <- ab2[1]
  bl <- abmix[2]
  b1 <- ab1[2]
  b2 <- ab2[2]
  P1 <- Pestimate(temperature,Tc1,Pc1,Wc1)
  P2 <- Pestimate(temperature,Tc2,Pc2,Wc2)
  P <- x * P1 + (1-x) * P2
  y1 <- 1 - round(x,1)
  y2 <- 1 - y1
  e <- 0.0001
   repeat
   {
     vl <- volsrk(temperature,P,al,bl,1)
     phil1 <- phisrk(temperature,P,vl,al,bl,a11,al12,b1,x)
     phil2 <- phisrk(temperature,P,vl,al,bl,al12,a22,b2,x)
     # Gas
     abmix <- calcabmix(a11,b1,a22,b2,k12,y1)
     ag <- abmix[1]
     ag12 <- abmix[3]
     bg <- abmix[2]
     vg <- volsrk(temperature,P,ag,bg,-1)
     phig1 <- phisrk(temperature,P,vg,ag,bg,a11,ag12,b1,y1)
     phig2 <- phisrk(temperature,P,vg,ag,bg,ag12,a22,b2,(1-y2))
     K1 <- phil1 / phig1
     K2 <- phil2 / phig2
     S <- x * K1 + (1-x)*K2
     if (abs(S-1) <= e){break}
     else
     {
      P <- P * S
      y1 <- (K1*x) / S
      y2 <- (K2*(1-x)) / S
     }
   }
  return(c(y1,P))#,phil1,phil2,phig1,phig2,vl,vg))
}

# calculate y, T, phi*, v* for a binary mixture
mixturecalct <- function(x,substance1,substance2,k12,P)
{
  substance <- lookupcritical(substance1)
  Tc1 <- as.numeric(substance[2])
  Pc1 <-  as.numeric(substance[3])
  Wc1 <-  as.numeric(substance[4])

  substance <- lookupcritical(substance2)
  Tc2 <- as.numeric(substance[2])
  Pc2 <-  as.numeric(substance[3])
  Wc2 <-  as.numeric(substance[4])

  temperature1 <- Testimate(P,Tc1,Pc1,Wc1)
  temperature2 <- Testimate(P,Tc2,Pc2,Wc2)
  temperature <- x * temperature1 + (1-x) * temperature2
  print(temperature)
  temperature<-314.6791
  Tr1 <- temperature / Tc1
  Tr2 <- temperature / Tc2
  ab1 <- calcab(Tc1,Pc1,Wc1,Tr1,temperature)
  ab2 <- calcab(Tc2,Pc2,Wc2,Tr2,temperature)
  # Liquid
  abmix <- calcabmix(ab1[1],ab1[2],ab2[1],ab2[2],k12,x)
  al <- abmix[1]
  a11 <- ab1[1]
  al12 <- abmix[3]
  a22 <- ab2[1]
  bl <- abmix[2]
  b1 <- ab1[2]
  b2 <- ab2[2]
  y1 <- 1 - round(x,1)
  y2 <- 1 - y1
  e <- 0.0001
 # repeat
#  {
    vl <- volsrk(temperature,P,al,bl,1)
    phil1 <- phisrk(temperature,P,vl,al,bl,a11,al12,b1,x)
    phil2 <- phisrk(temperature,P,vl,al,bl,al12,a22,b2,x)
    # Gas
    abmix <- calcabmix(a11,b1,a22,b2,k12,y1)
    ag <- abmix[1]
    ag12 <- abmix[3]
    bg <- abmix[2]
    vg <- volsrk(temperature,P,ag,bg,-1)
    phig1 <- phisrk(temperature,P,vg,ag,bg,a11,ag12,b1,y1)
    phig2 <- phisrk(temperature,P,vg,ag,bg,ag12,a22,b2,(1-y2))
    K1 <- phil1 / phig1
    K2 <- phil2 / phig2
    S <- x * K1 + (1-x)*K2
    print(S)
    if (abs(S-1) <= e){break}
   # else
  #  {
   #   temperature <- temperature + S
    #  y1 <- (K1*x) / S
    #  y2 <- (K2*(1-x)) / S
  #  }
#  }
  return(c(y1,temperature))#,phil1,phil2,phig1,phig2,vl,vg))
}
# end of functions
#
# global environment

# critical parameters for a couple of species
substances<-read.csv("~/Dropbox/MA-Thesis/Fuels/calculations/values/criticalvalues.csv",head=T,sep=",")
# universal gas constant
R <- 8.31433

# example 1: ethanol(17)-water(16)-mixture at 400 K over a range of liquid fractions (0 <= x <= 1)
# k12 <- -0.1# binary interaction parameter
# x <- seq(0,1,0.01) # x-sequence
# xyp <- c() # empty vector
# for (i in 1:length(x))
# {
#   xyp <- c(xyp,x[i],mixturecalc(x[i],17,16,k12,300)) # vector filled line per line with the results
# }
# srkXYP <- matrix(xyp,nrow=length(x),ncol=3,byrow=T) # matrix generated from the result vector
# plot(XYP[,1],XYP[,2],type='l',lty=4,lwd=2.5,col='orange1') # simple XY-line plot
# lines(XYP[,1],XYP[,1]) # a 45° line added to the plot
# # for a PXY plot: plot(XYP[,1],XYP[,3],type='l')
# #                 lines(XYP[,2],XYP[,3])
# rm(k12,x,xyp,i) # now useless values are removed from the global environment