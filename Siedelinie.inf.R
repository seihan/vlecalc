# load functions
source('compcheck.R')
source('Antoine.R')
source('Tboil.R')
Siedelinie.inf <- function(components,frac,P)
{
  #  Author:  Hannes Seidel  Hannes.Seidel@tu-cottbus.de   
  #        			                                          
  #  Name:	Siedelinie.inf		Version:	0.2 			                  
  #
  #  Estimate the bubble- and dewpoint for a mixture with given concentrations.                 
  #                                                         
  #  Txvle.inf (components,cofrac,u,r,q,system pressure)     
  #                                               
  #  components = [components,7] matrix of Antoine constants
  #  frac = mole fraction of each component
  #  u = [components,components] matrix                                       
  #  r,q = component dim vector                                     
  #                                                         
  #  returns the bubble- and dewpoint
  
  noc <- length(components) # count number of components
  if (length(frac) != noc){ # check: # substances == # fractions
    stop('Number of components differ their mole fractions. Abort.')
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (sum(frac)!=1){
    stop('Sum of mole fractions not equal to 1. Abort.')
  }

  # load Antoine constants
  complist <- c()
  for (i in 1:length(components)){
    complist <- rbind(complist,compcheck(components[i]))
  }
  # add mole fractions to complist
  complist <- cbind(complist,frac)
 
  # estimate boiling points for each component under ideal assumptions
  boil.T <- c()
  for (i in 1:length(components)){
    boil.T <- c(boil.T,Tboil(components[i],P))
  }
  # add boiling points to Complist
  complist <- cbind(complist,boil.T)
  # sort list by boiling points
  complist <- complist[order(complist[,9]), ]
  # 1. bubble point estimation
  e <- 0.00001 # best assumption, 0.00001 creates failures
  c1 <- 0 # iteration control
  Temperature <- complist[1,9] # start point (ideal)
  x <- complist[,8]
  repeat{
    #act <- uniquac(2,c(x[i],1-x[i]),u,r,q,Temperature) # activity coefficients
    Ps <- Antoine(complist[,2],complist[,3],complist[,4],Temperature)
    Pc <- sum(Ps * x) # calculated pressure
    y <- Ps * x / P
    if ((abs(sum(y)-1) < e) || (c1 > 999)){break}
    else if(sum(y) < 1){Temperature <- Temperature + sum(y)*0.1}
    else {Temperature <- Temperature - sum(y)*0.1}
    c1 <- c1 + 1 # count the iteration per x step
  }
  T0 <- Temperature # bubble point
  # 2. dewpoint estimation
  c2 <- 0
  repeat{
    y <- complist[,8] # vapour mole fraction at dewpoint
    Ps <- Antoine(complist[,2],complist[,3],complist[,4],Temperature) # list with saturated pressures
    x <-  y * P / Ps # liquid mole fraction at dewpoint
    Pc <- sum(Ps * x) # calculated pressure
    if ((abs(sum(x)-1) < e) || (c2 > 999)){break}
    else if(sum(x) < 1){Temperature <- Temperature - sum(x)*0.1}
    else {Temperature <- Temperature + sum(x)*0.1}
    c2 <- c2 + 1 # count the iteration per x step
  }
  T100 <- Temperature
  # 3. VLE calculation from T0 to T100 at P
  b <- y # original concentration -> solution vector
  L <- 1
  Llist <- c() # initiate list for liquid fraction
  c3list <- c() # list for iteration control
  Trange <- seq(T0,T100,1)
  for (i in 1:length(Trange)){
    c3 <- 0
    Ps <- Antoine(complist[,2],complist[,3],complist[,4],Trange[i]) # list with saturated pressures
    K <- Ps / P # Equilibrium constants
    repeat{
      xf <- L*(1-K) + K  # x factors for linear equation system
      A <- matrix(0,noc,noc) # initiate matrix
      for (i in 1:noc){
        A[i,i] <- xf[i] # fill the matrix
      }
      x <- solve(A,b) # solve the system
      if ((abs(sum(x) - 1) < e) && (abs(sum(K * x) - 1) < e) || (c3 > 999)){break}
      else if (abs(sum(K*x) < 1)){
        L <- L + e
      }
      else{
        L <- L - 0.001
      }
      c3 <- c3 + 1
    }
    Llist <- rbind(Llist,L)
    c3list <- rbind(c3list,c3)
  }
  colnames(c3list) <- 'iterations'
  Tl <- cbind(Trange,Llist*100)
  Tl <- rbind(Tl,c(T100,0))
  colnames(Tl) <- c('T [K]','% liquid')
  T0100 <- c(T0,T100) # T.bubble, T.dew
  iT0100 <- c(c1,c2) # iterations
  return(list(T0100,iT0100,c3list,Tl))
}
## Examples, uncomment for using
#Siedelinie.inf(c(5,1,30),c(0.2,0.4,0.4),1.013e+5)