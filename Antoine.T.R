Antoine.T<-function(A,B,C,P){
  # Calculates the saturation temperature with the Antoine equation
  # Ts = B / (A - log10(P)) - C
  Ts = B / (A - log10(P*1e-5)) - C;
  return(as.numeric(Ts))
}