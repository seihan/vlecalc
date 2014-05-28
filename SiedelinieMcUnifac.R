# Author:  Hannes Seidel  Hannes.Seidel@tu-cottbus.de   
SiedelinieMcUnifac <- function(components,frac,unu,aij,kP){
  # Estimate the bubble- and the dewpoint for a mixture with the UNIFAC
  # model and calculate the vapour liquid equilibrium between under ideal
  # assumptions.
  #
  # Args:                                             
  # components: numbers for components
  # frac: mole fraction of each component
  # unu: UNIFAC functional groups matrix                                       
  # aij: UNIFAC interaction parameters matrix                                     
  # kP: System pressure                                                        
  #
  # Returns:
  # The bubble- and dewpoint and the vle data between
  # load functions
  source('Compcheck.R')
  source('AntoineP.R')
  source('AntoineT.R')
  source('Tboil.R')
  source('Unifac.R')
  source('printf.R')
  Percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  kR <- 8.314
  noc <- length(components) # count number of components
  if (length(frac) != noc){ # check: # substances == # fractions
    stop('Number of components differ their mole fractions. Abort.')
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (sum(frac)!=1){
    stop('Sum of mole fractions not equal to 1. Abort.')
  }
  # add mole fractions to complist
  complist <- cbind(Compcheck(components), frac)
  Ps <- function(temperature, ...) {
    Antoine.P(complist[, 2], complist[, 3], complist[, 4], temperature)
  } # list with saturated pressures
  # 1. bubble point estimation
  e <- 0.0001  # accuracy
  c1 <- 0  # iteration control
  x <- complist[, 9]
  temperature <- sum(x * Antoine.T(complist[, 2],
                                   complist[, 3],
                                   complist[, 4],
                                   kP))
  rat <- 0
  y <- c()
  while (c1 < 100) {
    act <- Unifac(x, unu, aij, temperature)
    pCalc <- sum(act * Ps(temperature) * x)  # calculated pressure
    rat <- kP / pCalc
    if (abs(rat - 1) < e) 
      break
    y <-  Ps(temperature) * x / pCalc
    v <- kR * temperature / pCalc
    f <- act * Ps(temperature) * exp((v * x) * (pCalc - kP)/(kR * temperature))
    psNew <- (y * act * kP/ (x * f)) * Ps(temperature)
    temperature <- sum(x*Antoine.T(complist[, 2],
                                   complist[, 3],
                                   complist[, 4],
                                   psNew))
    c1 <- c1 + 1  # count the iteration
  }
  pCalc1 <- pCalc
  t0 <- temperature  # bubble point
  printf('The bubble point temperature is %3.3f K (%3.3f°C).\n',
         t0,
         t0-273.15)
  # 2. dewpoint estimation
  temperature <- sum(x * Antoine.T(complist[, 2],
                                   complist[, 3],
                                   complist[, 4],
                                   kP))
  act <- rep(1, length(y))
  tNew <- temperature + 50
  dT <- 50  # temperature difference
  c2 <- 0
  y <- complist[, 9] # vapour mole fraction at dewpoint
  objF <- function(temperature, kP, act, y){  # objective function
    kP - (1 / sum(y / (act * Ps(Temperature))))
  }
  while (c2 < 10) {
    tFormer <- temperature
    k <- kP / (act * Ps(temperature)) # equilibrium constant
    pCalc <- 1 / sum(y / (act * Ps(temperature)))
    temperature<-as.numeric((optimize(f=objF,  # minimize the objective function
                                      lower=tFormer,
                                      upper=tFormer+dT,
                                      kP=kP,
                                      act=act,
                                      y=y,
                                      tol=100))[1])
    tNew <- temperature
    dT <- tNew - tFormer
    x <- y / k
    act <- Unifac(x, unu, aij, temperature)
    rat <- pCalc / kP
    if (abs(rat - 1) < e){break}
    S <- sum(x)
    c2 <- c2 + 1 # count the iteration per x step
  }
  temp <- temperature
  c2 <- 0
  dT <- 10
  X <- function(...)(y * pCalc / Ps(temperature))
  objF <- function(temperature,...) {
    1 - sum(X(y, pCalc, Ps(temperature)))
  }
    while(c2 < 10){
    tFormer <- temperature
    pCalc <- 1 / sum(y / (act * Ps(temperature)))
    temperature <- as.numeric(optimize(objF,  # minimize the objective function
                                       lower=tFormer,
                                       upper=tFormer + dT,
                                       temperature=tFormer,
                                       kP=pCalc,
                                       Ps=Ps(Temperature),
                                       y=y,
                                       tol=1000)[1])
    tNew <- temperature
    dT <- tNew - tFormer
    act <- Unifac(X(y, pCalc, Ps(Temperature)), unu, aij, temperature)
    c2 <- c2 + 1
  }
  t100 <- (temp + temperature) / 2
  if(t100 > 1000){
    printf('The dew point temperature is %3.3f K (%3.3f°C). This is to much. Abort.',
           t100,
           t100-273.15)
    stop('Temperature exceed vality range.')
  }
  printf('The dew point temperature is %3.3f K (%3.3f°C).\n\nStarting VLE calculation ...\n',
         t100,
         t100-273.15)
  # 3. VLE calculation from T0 to T100 at P
  b <- y # original concentration -> solution vector
  l <- 1 # liquid part
  lList <- c() # initiate list for liquid fraction
  c3List <- c() # list for iteration control
  tRange <- seq(T0,T100,1)
  #e <- 0.01
  for (i in 1:length(tRange)){
    printf('%s ', Percent(i / length(Trange)))
    c3 <- 0
    #Ps <- Antoine(complist[,2],complist[,3],complist[,4],Trange[i]) # list with saturated pressures
    k <- Ps(tRange[i]) / kP # Equilibrium constants
    repeat {
      xf <- l * (1 - k) + k  # x factors for linear equation system
      A <- matrix(0, noc, noc)  # initiate matrix
      for (i in 1:noc) {
        A[i, i] <- xf[i] # fill the matrix
      }
      x <- solve(A, b)  # solve the system
      if ((abs(sum(x) - 1) < e) && (abs(sum(k * x) - 1) < e) || (c3 > 999)) {
        break
      } else if (abs(sum(K*x) < 1)) {
        l <- l + e
      }
      else {
        l <- l - 0.001
      }
      c3 <- c3 + 1
    }
    lList <- rbind(lList, l)
    c3List <- rbind(c3list, c3)
  }
  colnames(c3List) <- 'iterations'
  temLiquid <- cbind(tRange, lList*100)
  temLiquid <- rbind(temLiquid,c(t100,0))
  colnames(temLiquid) <- c('T [K]','% liquid')
  t0100 <- c(t0,t100) # T.bubble, T.dew
  iT0100 <- c(c1,c2) # iterations
  return(list(T0100, c(pCalc1, pCalc), iT0100, c3List, temLiquid))
}
## Examples, uncomment for using2
#Siedelinie.inf(c(5,1,30),c(0.2,0.4,0.4),1.013e+5)
#Siedelinie.mc.real(c(13,15,14,16,17,24,23,27,30),c(0.099,0.091,0.198,0.08,0.01,0.314,0.06,0.05,0.098),1,2,1.013e+5)
#n-heptylbenzene,naphthalene,tetralin,n-hexadecane,cyclohexane,n-octadecane,phenanthrene,n-tetradecane,n-heneicosane
#udat<-unifactool2()
#unu <- udat[[1]]
#aij <- udat[[2]]
FACE6<-Siedelinie.mc.real(c(30,24,23,14,17,15,27,13,16),c(0.098,0.314,0.060,0.198,0.010,0.091,0.050,0.099,0.080),unu,aij,1.013e+5)
#   boil.T <- c()
#   for (i in 1:noc){
#     boil.T <- c(boil.T,Tboil(components[i],P))
#   }
#   # add boiling points to Complist
#   complist <- cbind(complist,boil.T)
#   # sort list by boiling points
#   complist <- complist[order(complist[,10]), ]
#   print(complist[,9])
#   return(complist)
pdf("vapline-FACE6.pdf")
plot(100-FACE6[[5]][,2],FACE6[[5]][,1],type='l',lwd=2,main='Vaporization Line - FACE#6',xlab='Evaporated fraction [%]',ylab='T [K]')
#legend(0.05,1,c("Experiment (DDB)","Ideal","Van Laar","Margules","SRK (k12 = -0.1)","UNIQUAC (PHI = 1)"),lty=c(NA,1,6,4,2,2),lwd=c(1.5,2,2,2,2),pch=c(8,NA,NA,NA,NA,NA),col=c("black","gray66","blue2","red","firebrick","seagreen"),cex=0.7,merge=FALSE)
dev.off()

