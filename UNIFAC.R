UNIFAC = function(fractions=c(),
                  unu=c(), 
                  aij=c(),
                  temperature=NULL){
  # function to calculate the activity coefficients 
  # with UNIFAC
  # arguments:
  # x: vector with molar fractions
  # unu: UNIFAC surface parameters and component matrix
  # aij: UNIFAC interaction parameters matrix [k,k]
  # temperature: [K]
  #
  # return:
  # activity coefficients
  nos = length(fractions);# count number of substances
  e = 1e-7;
  # check for mole fraction consistent (sum(xi = 1)?)
  if (abs(sum(fractions) - 1) > e){
    stop('Sum of mole fractions not equal to 1. Abort.')
  }
  x = fractions;
  # 1. create UNIFAC matrices
  r = rep(0, nos);  # r volume contribution, c(r[1], ... ,r[n])
  for (i in 1:nos){
    r[i] = sum(unu[, 3] * unu[, 4 + i]);
  }
  q = rep(0, nos);  # q surface fractions, c(q[1], ... , q[n])
  for (i in 1:nos){
    q[i] = sum(unu[, 4]*unu[, 4 + i]);
  }
  nog = nrow(unu); # number of groups
  eij = matrix(0, nog, nos); # surface fractions for
  for (i in 1:nos){ # each subgroup in each molecule
    for (j in 1:nog){
      eij[j, i] = unu[j, 4 + i] * unu[j, 4] / q[i];
    }
  }
  # 1. summation
  # determine the overall surface area fraction 
  # for each subgroup in the mixture
  S = rep(0, nog);
  for (i in 1:nog){
    S[i] = sum(x * q * eij[i, ]) / sum(x * q);
  }
  tau = exp(-aij / temperature);
  beta = matrix(0, nos, nog);
  for (i in 1:nos){
    for (j in 1:nog){
      beta[i, j] = sum(eij[, i] * tau[, j]);
    }
  }
  s = rep(0, nog);
  for (i in 1:nog){
    s[i] = sum(S * tau[, i]);
  }
  J = rep(0, nos)
  for (i in 1:nos){
    J[i] = r[i] / sum(x * r);
  }
  L = rep(0, nos);
  for (i in 1:nos){
    L[i] = q[i] / sum(x * q);
  } # 2. combinatory part
  actc = rep(0, nos);
  for (i in 1:nos){
    actc[i] = 1 - J[i] + log(J[i]) - 5 * q[i] * 
      (1 - J[i]/ L[i] + log(J[i] / L[i]));
  } # 3. residual part
  actr = rep(0,nos);
  for (i in 1:nos){
    actr[i] = q[i] * (1- sum(S * beta[i,] / s -
                               eij[,i] * log(beta[i,] / s)));
  } # 4. sum the parts
  act = exp((actc + actr));
  return(act);
}