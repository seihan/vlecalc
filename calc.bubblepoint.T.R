# Author:  Hannes Seidel  Hannes.Seidel@b-tu.de
calc.bubblepoint.T <- function(Substances, Fractions, Pressure){
  # function to calculate the bubblepoint of multicomponent mixtures
  # 1. Bubble Point Estimation
  # Arguments:
  #   substances: Vector with names of substances c('ethanol','n-decan')
  #   fractions: a vector with mole fraction of each substance
  #   P: System pressure
  #
  # Returns:
  #   bubblepoint temperature [K],
  #   activity coefficients,
  #   gas fractions (y)
  #   iteration count (c)
  #
  # Load functions:
  source('Antoine.P.R')
  source('Antoine.T.R')
  source('sub.check.R')
  source('UNIFAC.gen.R')
  source('UNIFAC.R')
  Ps <- function(Temperature, ...) {
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature)
  }    
  kR <- 8.314 # universal gas constant
  P = Pressure
  Uvalues <- UNIFAC.gen(Substances) # load or generate UNIFAC values
  unu = Uvalues[[1]]
  aij = Uvalues[[2]]
  nos <- length(Substances) # count number of components
  if (length(Fractions) != nos){ # check: # substances == # fractions
    stop('Number of components differ their mole fractions. Abort.')
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(sum(Fractions) - 1) > 0.001 ){
    print(Fractions)
    stop('Sum of mole fractions not equal to 1. Abort.')
  }
  sublist <- NULL
  for (i in 1:length(Substances)){
    tmp <- sub.check(Substances[i])
    sublist <- rbind(sublist, tmp[1,])
  }
  # add mole fractions to sublist
  sublist <- cbind(sublist, Fractions)
  e <- 1e-7  # accuracy
  c <- 0  # iteration count
  x <- Fractions
  Temperature <- sum(x * Antoine.T(sublist[, 2],
                                   sublist[, 3],
                                   sublist[, 4],
                                   Pressure))
  rat <- 0
  y <- c()
  while (c < 1000) { # emergency break
    act <- UNIFAC(x, Uvalues[[1]], Uvalues[[2]], Temperature)
    if(is.na(act[1])) break
    Pcalc <- sum(act * Ps(Temperature) * x)  # calculated pressure
    rat <- Pressure / Pcalc
    if (abs(rat - 1) < e) break
    y <-  Ps(Temperature) * x / Pcalc
    v <- kR * Temperature / Pcalc
    f <- act * Ps(Temperature) * exp((v * x) * (Pcalc - Pressure)/(kR * Temperature))
    PsNew <- (y * act * P / (x * f)) * Ps(Temperature)
    Temperature <- sum(x * Antoine.T(sublist[, 2],
                                     sublist[, 3],
                                     sublist[, 4],
                                     PsNew))
    c <- c + 1  # count the iteration
  }
#   printf('\nThe bubble point temperature is %3.3f K (%3.3f°C).\n',
#          temperature,
#          temperature-273.15)
  Temperature
#  return(list(temperature,act,y,c))
}