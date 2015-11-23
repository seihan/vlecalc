calc.boiling.T <- function (Substance, P){
  Substance <- sub.check(Substance)
  temperature <- Antoine.T(Substance[1,2], Substance[1,3], Substance[1,4], P)
  return(temperature)
}