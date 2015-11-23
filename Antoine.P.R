Antoine.P = function(A,
                     B,
                     C,
                     Temperature){
  # Calculates the saturation pressure with the Antoine equation
  # p^s = 10^(A - B / (C + T[K])
  Ps = (10^(A - B / (C + Temperature)))*10^5; # constances in [bar]
  return(as.numeric(Ps))
}