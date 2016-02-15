os.boilingline = function(substances=NULL,
                          fractions=NULL,
                          pressure=NULL){
  source('calc.bubblepoint.T.R');
  source('Antoine.P.R');
  source('Antoine.T.R');
  source('SRK.R');
  Percent = function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  Ps = function(Temperature, ...) {# list with saturated pressures
    Antoine.P(sublist[, 2], sublist[, 3], sublist[, 4], Temperature)
  }
  printf = function(...) cat(sprintf(...));
  e <- 1e-7  # accuracy
  nos = length(substances); # count number of substances
  if(is.null(pressure)){
    stop('No pressure [Pa] specified. Abort.');
  }
  if (nos == 0){
    stop('No substance specified. Abort.');
  }
  if (length(fractions) != nos){ # check: N substances == N fractions
    stop('Number of components differ their mole fractions. Abort.');
  }
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(1 - sum(fractions)) > e){
    stop('Sum of mole fractions not equal to 1. Abort.');
  }
  Substances = rep(list(''), nos); # initiate the list of lists
  for(i in 1:nos){
    Substances[[i]] = Substance(substances[i]); # get properties
    Substances[[i]]$Fraction = fractions[i]; # add fractions
    A = Substances[[i]][[1]]$A[1];
    B = Substances[[i]][[1]]$B[1];
    C = Substances[[i]][[1]]$C[1];
    Tsat = Antoine.T(A, B, C, pressure) # calc boiling temperature
    Substances[[i]]$Tsat = Tsat;
  }
  names(Substances) = substances;
  oSubstances = Substances[order(-sapply(Substances, function(tsat) tsat$Tsat))]; # sort by Tsat
  M = sum(sapply(Substances, function(molarmass) molarmass$MolarMass)); # sum molar mass
  stop("Under construction");
}