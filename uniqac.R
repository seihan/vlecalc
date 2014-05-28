##########################################################
##                      UNIQUAC                         ##
##              Universal Quasichemical                 ##
##   -----------------------------------------------    ##
##                                                      ##
##   -----------------------------------------------    ##
##                                                      ##
##   input:  nkomp,x[i],a[i,j],r[i],q[i],T[K]           ##
##           T                                          ##
##   output: activity coefficients[i]                   ##
##                                                      ##
##########################################################

uniqac <- function (nkomp,x,aij,r,q,temperature)
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