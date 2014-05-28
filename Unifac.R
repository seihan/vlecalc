#  Author:  Hannes Seidel  Hannes.Seidel@tu-cottbus.de   
Unifac <- function(x, unu, aij, temperature) {
  # function to calculate the activity coefficients using UNIFAC
  # 
  # Args:
  # x: vector with molar fractions, length(x) == noc !
  # unu: UNIFAC surface parameters and component matrix
  # aij: UNIFAC interaction parameters matrix [k,k]
  # temperature: Temperature [K]
  #
  # Returns:
  # Activity coefficients
  noc <- length(x)  # number of components
  # 1. create UNIFAC matrices
  r <- rep(0, noc)  # r volume contribution, c(r[1], ... ,r[n])
  for (i in 1:noc) {
    r[i] <- sum(unu[, 3] * unu[, 4 + i])
  }
  q <- rep(0, noc)  # q surface fractions, c(q[1], ... , q[n])
  for (i in 1:noc) {
    q[i] <- sum(unu[, 4]*unu[, 4 + i])
  }
  nog <- nrow(unu) # number of groups
  eij <- matrix(0, nog, noc) # surface fractions for each subgroup in each molecule
  for (i in 1:noc){
    for (j in 1:nog){
      eij[j, i] <- unu[j, 4 + i] * unu[j, 4] / q[i]
    }
  }
  # 1. summation
  # determine the overall surface area fraction for each subgroup in the mixture
  S <- rep(0, nog)
  for (i in 1:nog) {
    S[i] <- sum(x * q * eij[i, ]) / sum(x * q)
  }
  tau <- exp(-aij / temperature)
  beta <- matrix(0, noc, nog)
  for (i in 1:noc) {
    for (j in 1:nog) {
      beta[i, j] <- sum(eij[, i] * tau[, j])
    }
  }
  s <- rep(0, nog)
  for (i in 1:nog) {
    s[i] <- sum(S * tau[, i])
  }
  J <- rep(0, noc)
  for (i in 1:noc) {
    J[i] <- r[i] / sum(x * r)
  }
  L <- rep(0, noc)
  for (i in 1:noc) {
    L[i] <- q[i] / sum(x * q)
  }
  # 2. combinatory part
  actc <- rep(0, noc)
  for (i in 1:noc) {
    actc[i] <- 1 - J[i] + log(J[i]) - 5 * q[i] * (1 - J[i]/L[i] + log(J[i]/L[i]))
  }
  # 3. residual part
  actr <- rep(0,noc)
  for (i in 1:noc){
    actr[i] <- q[i] * (1- sum(S * beta[i,] / s - eij[,i] * log(beta[i,] / s)))
  }
  # 4. sum the parts
  act <- exp((actc + actr))
  return(act)
}
# Example
# ulist <- unifactool() # interactive user input
# unu <- ulist[[1]] # UNIFAC surface parameters and component matrix
# aij <- ulist[[2]] # UNIFAC interaction parameters matrix [k,k]
# unifac(c(0.4,0.6),unu,aij,308.15) # returns the activity coefficients for two substances