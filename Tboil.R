source('Compcheck.R')
Tboil<-function(comp,P){
  # Estimate the boiling temperature under ideal conditions, using the Antoine equation
  comp<-Compcheck(comp)
  Temperature <- 300
  e <- 0.0001
  c <- 0 # iteration control and emergency break
  repeat
  {    
    Ps<-Antoine(comp[,2],comp[,3],comp[,4],Temperature)
    rat <- Ps / P
    if ((abs(1 - rat) < e) || (c > 9999)){break}
    else if (rat > 1){Temperature <- Temperature - rat}
    else if (rat < 1e-4){Temperature <- Temperature + rat * 1e+6} # in both directions and take care to stay in range
    else {Temperature <- Temperature + rat}
    c <- c + 1
  }
  return(c(Temperature))
}