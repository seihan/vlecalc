Siedlinie.inf <- function (components,u,r,q,P) # main function, co functions are below
{
  #  Author:  Hannes Seidel  Hannes.Seidel@tu-cottbus.de   
  #        			                                          
  #  Name:	Txvle.inf		Version:	0.2 			                  
  #
  #  This function calculate the molar fractions of the     
  #  liquid and gas phase for binary mixtures.              
  #  It starts at the boiling point of the highboiler.       
  #  The saturated pressure and the activity coefficient    
  #  are calculated using the Antoine equation and the      
  #  UNIQUAC algorithm.                                     
  #  The gas phase is ideal.                                
  #                                                         
  #  Txvle.inf (components,u,r,q,system pressure)     
  #                                               
  #  components = [components,7] matrix of Antoine constants
  #  u = [components,components] matrix                                       
  #  r,q = component dim vector                                     
  #                                                         
  #  returns a resultant matrix with T,xi,yi,acti,setps
  
  ## x -> psat, act @ T.boiling -> sum (x * psat * act) == P ? -> T(P,x), y(P,T,x)
  
  # load Antoine constants
  complist <- c()
  for (i in 1:length(components)){
    complist <- c(complist,compcheck[i]))
  }
  Complist <- matrix(complist,nrow=length(components),ncol=7,byrow=T)
  # estimate boiling points for each component under ideal assumptions
  tboil <- c()
  for (i in 1:length(components)){
    tboil <- c(tboil,Tboil(components[i],P))
  }
  # add boiling points to Complist
  Complist <- cbind(Complist,tboil)
  # sort list by boiling points
  Complist <- gnomesort.matrix(Complist,8)
  e <- 0.001 # best assumption, 0.0001 creates failures
  x<-seq(0,1,0.01) # x vector
  txy<-c() # initiate the resultant vector
  Temperature <- T0 # start point (ideal)
  for (i in 1:length(x))
  {
    c <- 0 # iteration control (emergency break)
    repeat
    {
      Ps1 <- Antoine(comp1[2],comp1[3],comp1[4],Temperature) # saturated pressure
      Ps2 <- Antoine(comp2[2],comp2[3],comp2[4],Temperature)
      act <- uniquac(2,c(x[i],1-x[i]),u,r,q,Temperature) # activity coefficients
      rat <- (x[i]*Ps1*act[1] + (1-x[i])*Ps2*act[2])/ P # Pcalc / P - ratio
      if ((abs(rat - 1) < e)|| (c > 999)) {break} # solution criteria and emergency break
      else if (rat < 1) {Temperature <- Temperature + rat} # decrease temperature with the rat
      else {Temperature <- Temperature - rat} # increase
      c <- c + 1 # count the iteration per x step
    }
    y <- x[i]*Ps1*act[1]/(x[i]*Ps1*act[1] + (1-x[i])*Ps2*act[2]) # gas fraction y1 <- x1*Ps1*act1 / P
    txy <- c(txy,Temperature,x[i],y,act,c) # fill the resultant vector
  }
  TXY<-matrix(txy,nrow=length(x),ncol=6,byrow=T) # transform to resultant matrix
  colnames(TXY)<-c("T [K]","x(T,P)","y(T,P,x)","act 1","act 2","steps") # set the column names
  return(TXY)
}

uniquac <- function (nkomp,x,aij,r,q,temperature)
{
  # calculating the combinatorially part
  Sri <- 0
  Sqi <- 0
  qx <- rep(NA,nkomp)  
  lngc <- rep(NA,nkomp)
  for (i in 1:nkomp)
  {
    Sri <- Sri + r[i] * x[i]
    qx[i] <- q[i] * x[i]
    Sqi <- Sqi + qx[i]
  }
  for (i in 1:nkomp)
  {
    rr <- r[i] / Sri
    qq <- q[i] / Sqi
    rqi <- rr / qq
    lngc[i] <- 1 - rr + log(rr) - 5 * q[i] * (1 - rqi + log(rqi))
  }
  
  # calculating the rest part
  tau <- matrix(0, nkomp, nkomp)
  lngr <- rep(NA,nkomp)
  gam <- rep(NA,nkomp)
  for (i in 1:nkomp)
  {
    for (j in 1:nkomp)
    {
      tau[i,j] <- exp(-aij[i,j] / temperature)
    }
  }
  for (i in 1:nkomp)
  {
    sk1 <- 0
    sk2 <- 0
    for (j in 1:nkomp)
    {
      sk1 <- sk1 + qx[j] * tau[j,i]
      sk3 <- 0
      for (k in 1:nkomp)
      {
        sk3 <- sk3 + qx[k] * tau[k,j]
      }
      sk2 <- sk2 + (qx[j] * tau[i,j]) / sk3
    }
    lngr[i] <- (1 - log(sk1 / Sqi) - sk2) * q[i]
    gam[i] <- exp(lngc[i] + lngr[i])
  }
  return(gam)
}

compcheck<-function(comp)
{
  for (i in 1:nrow(compvals))
  {
    if (comp==i)
    {
      comp<-compvals[i,]
      break
    }
  }
  return(comp)
}

# Estimate the boiling temperature
Tboil<-function(comp,P)
{
  comp<-compcheck(comp)
  Temperature <- 300
  e <- 0.0001
  i <- 1
  repeat
  {    
    Ps<-Antoine(comp[,2],comp[,3],comp[,4],Temperature)
    rat <- P / Ps
    if ((rat - 1 < e)|| (i == 100)){break}
    else if (rat < 1){Temperature <- Temperature - rat}
    else if (rat > 1e+2){Temperature <- Temperature + rat/1e+4} # in both directions and take care to stay in range
    else {Temperature <- Temperature + rat}
    i <- i + 1
  }
  return(Temperature)
}
## ANTOINE    ##
# Calculates the saturation pressure with the Antoine equation
# p^s = 10^(A - B / (C + T[K])
Antoine<-function(A, B, C, Temperature)
{
  psat <- (10^(A - B / (C + Temperature)))*10^5
  return(as.numeric(psat))
}
## gnomesort
# Sort vector
gnomesort <- function(vec){
  i <- 1  
  while (i <= length(vec)){
    if (i == 1 || vec[i-1] <= vec[i]){
      i <- i + 1
    }
    else{
      tmp <- vec[i]
      vec[i] <- vec[i-1]
      i <- i - 1
      vec[i] <- tmp
    }
  }
  return(vec)
}
## gnomesort.matrix
# Sort matrix by specific column
gnomesort.matrix <- function(mat,col){
  i <- 1  
  while (i <= nrow(mat)){
    if (i == 1 || mat[i-1,col] <= mat[i,col]){
      i <- i + 1
    }
    else{
      tmp <- mat[i,]
      mat[i,] <- mat[i-1,]
      i <- i - 1
      mat[i,] <- tmp
    }
  }
  return(mat)
}

# Antoine values
compvals<-read.csv("~/Dropbox/MA-Thesis/Fuels/calculations/values/comvals.csv",head=T,sep=",")

## Examples, uncomment for using

# Example EtOH-H2O, # EtOH: compvals[4,], H20: compvals[9,], P = Patm = 101325 Pa

# u <- matrix(c(0,175.47,-21.42,0),2)
# r <- c(2.1055,0.92)
# q <- c(1.972,1.4)
# P <- 101325
TXY <- Txvle(4,9,u,r,q,P*2)
