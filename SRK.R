SRK = function(pressure=NULL,
               temperature=NULL,
               x=NULL,
               y=NULL,
               Tc=NULL, # T critical [K]
               Pc=NULL, # P critical [bar]
               Ac=NULL, # acentric factor
               wij=NULL){# binary interaction parameter
  P = pressure;
  R = 0.0831433;
  PI = 3.141592653589793;
  e = 1e-7; # accuracy
  ABSRK = function(temperature=NULL,
                   x=NULL, 
                   nos=NULL, ...){
    if(is.null(temperature))stop('temperature = NULL')
    if(is.null(x))stop('x = NULL')
    am = alpha = bi = rep(0, nos);
    aij = matrix(NA, nos, nos);
    a = b = 0;
    for (i in 1:nos){
      am[i] = 0.48 + 1.574 * Ac[i] - 0.176 * Ac[i]^2;
      alpha[i] = (1 + am[i] * (1 - sqrt(temperature / Tc[i])))^2;
      aij[i,i] = 0.42748 * R^2 * Tc[i]^2 / Pc[i] * alpha[i];
      a = a + (x[i] * x[i] * aij[i,i]);
      bi[i] = 0.08664 * R * Tc[i] / Pc[i];
      b = b + (x[i] * bi[i]);
    }
    if(nos == 1){
      result = list(am = am,
                    alpha = alpha,
                    aij = aij,
                    a = a,
                    bi = bi,
                    b = b);
    }
    else{
      if(is.null(wij)) wij = matrix(0.01, nos, nos); # binary interaction parameter k12
      for(i in 1:(nos-1)){
        for(j in (i+1):nos){
          aij[i,j] = sqrt(aij[i,i] * aij[j,j]) * (1 - wij[i,j]);
          aij[j,i] = aij[i,j];
          a = a + 2 * x[i] * x[j] * aij[i,j];
        }
      }
      result = list(am = am,
                    alpha = alpha,
                    aij = aij,
                    a = a,
                    bi = bi,
                    b = b);
    }
    return(result)
  } # end ABSRK
  VOLSRK = function(phase=NULL,
                    temperature=NULL,
                    P=NULL,
                    a=NULL,
                    b=NULL, ...){
    pstr = -(P * b^2 + R * temperature * b - a) / 3 / P - (R * temperature / 3 / P)^2;
    qstr = -(R * temperature / 3 / P)^3 - R * temperature * (P * b^2 + R * temperature * b - a) / 
           6 / P^2 - a * b / 2 / P;
    diskr = qstr^2 + pstr^3;
    if (diskr < 0){
      rstr = sign(qstr) * sqrt(abs(pstr));
      cosphi = qstr / rstr^3;
      phi = acos(cosphi);
      x1 = -2 * rstr * cos(phi / 3);
      x2 = 2 * rstr * cos((PI - phi) / 3);
      x3 = 2 * rstr * cos((PI + phi) / 3);
      if (phase == -1){
        v = max(x1,x2,x3);
        }
      else {
        v = min(x1,x2,x3);
      }
    }
    else{
      h1 = -qstr + sqrt(diskr);
      h2 = -qstr - sqrt(diskr);
      h3 = sign(h1);
      h4 = sign(h2);
      v = h3 * abs(h1)^(1/3) + h4 * abs(h2)^(1/3);
    }
    volsrk = v + R * temperature / 3 / P;
    return(volsrk)
  } # end VOLSRK
  PHISRK = function(temperature=NULL,
                    P=NULL,
                    v=NULL,
                    x=NULL,
                    a=NULL,
                    b=NULL,
                    bi=NULL,
                    aij=NULL,
                    nos=NULL, ...){
    h1 = log(v / (v - b));
    h2 = b / (v + b);
    h3 = log((v + b) / v);
    h4 = 1 / (R * temperature * b);
    h5 = a * h4 * (h3 - h2) / b;
    h6 = log(P * v / (R * temperature));
    phi = rep(0, nos);
    for (i in 1:nos){
      phi[i] = h1 + bi[i] / (v - b) - h6 + bi[i] * h5;
      S = 0;
      for (j in 1:nos){
        S = S + x[j] * aij[i,j];
      }
      phi[i] = phi[i] - 2 * S * h4 * h3;
      phi[i] = exp(phi[i]);
    }
    return(phi);
  } # end PHISRK
  Pestimate = function(temperature,Tc,Pc,Ac)
  {
    Temp = log(Pc);
    Temp1 = log(10) * (1 - Tc / temperature) * (7 + 7 * Ac) / 3;
    P = exp(Temp + Temp1);
    return(P);
  }
  Testimate = function(P,Tc,Pc,Ac)
  {
    Temp = log(P / Pc);
    Temp1 = 1 - Temp * 3 / (log(10) * (7 + 7 * Ac));
    temperature = Tc / Temp1;
    return(abs(temperature));
  }
  Kestimate = function(temperature,P,Tc,Pc,Ac)
  {
    Temp = log(Pc / P);
    Temp1 = log(10) * (1 - Tc / temperature) * (7 + 7 * Ac) / 3;
    K = exp(Temp + Temp1);
    return(K);
  }
  calc.bubblepoint.P = function(temperature, x, ...){
    P = Pestimate(temperature, Tc, Pc, Ac);
    P = sum(x * P);
    y = Kestimate(temperature, P, Tc, Pc, Ac) * x;
    y = y / sum(y)
    nos = length(x); # number of substances
    c = 0; # iteration control
    # 1. calc. a, b liquid phase
    abl = ABSRK(temperature,
                x,
                nos=nos);
    repeat{
      vl = VOLSRK(phase=1,    # 2. calc. v, phi liquid @ P
                  temperature,
                  P,
                  abl$a,
                  abl$b);
      phil = PHISRK(temperature,
                    P,
                    vl,
                    x,
                    abl$a,
                    abl$b,
                    abl$bi,
                    abl$aij,
                    nos=nos);
      abg = ABSRK(temperature, # 3. calc. a, b gas phase
                  y,
                  nos);
      vg = VOLSRK(phase=-1,  # 4. calc. v, phi gas @ P
                  temperature,
                  P,
                  abg$a,
                  abg$b);
      phig = PHISRK(temperature, 
                    P,
                    vg,
                    y,
                    abg$a,
                    abg$b,
                    abg$bi,
                    abg$aij,
                    nos=nos); 
      K = phil / phig;  # 5. calc. K = phi_l / phi_g, S = sum(x * K)
      S = sum(x * K);
      if ((abs(S - 1) < e) || (c > 99)){break;}
      else{
        c = c + 1;
        P = S * P;
        y = (K * x) / S;
        y = y / sum(y);
      }
    }
    result = list(pressure=P*1e+5,
                  y=y,
                  phil=phil,
                  phig=phig,
                  S=S,
                  iterations=c);
    return(result);
  } # end calc.bubblepoint.P
  calc.bubblepoint.T = function(P, x, ...){
    P = P * 1e-5;
    temperature = Testimate(P, Tc, Pc, Ac);
    temperature = sum(x * temperature);
    y = Kestimate(temperature, P, Tc, Pc, Ac) * x;
    y = y / sum(y)
    nos = length(x); # number of substances
    c = 0;
    repeat{
      abl = ABSRK(temperature, 
                  x,
                  nos);# 1. calc. a, b liquid phase
      vl = VOLSRK(phase=1,    # 2. calc. v, phi liquid @ P
                  temperature,
                  P,
                  abl$a,
                  abl$b);
      phil = PHISRK(temperature,
                    P,
                    vl,
                    x,
                    abl$a,
                    abl$b,
                    abl$bi,
                    abl$aij,
                    nos=nos);
      abg = ABSRK(temperature, # 3. calc. a, b gas phase
                  y,
                  nos);
      vg = VOLSRK(phase=-1,  # 4. calc. v, phi gas @ P
                  temperature,
                  P,
                  abg$a,
                  abg$b);
      phig = PHISRK(temperature, 
                    P,
                    vg,
                    y,
                    abg$a,
                    abg$b,
                    abg$bi,
                    abg$aij,
                    nos=nos); 
      K = phil / phig;  # 5. calc. K = phi_l / phi_g, S = sum(x * K)
      S = sum(x * K);
      if ((abs(S - 1) < e) || (c > 99)){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (1 - S) / S;
        temperature = temperature + temp;
        y = (K * x) / S;
        y = y / sum(y);
      }
    }
    result = list(temperature=temperature,
                  y=y,
                  phil=phil,
                  phig=phig,
                  S=S,
                  iterations=c);
    return(result);
  } # end calc.bubblepoint.T
  calc.dewpoint.P = function(temperature, y, ...){ # calculate dewpoint pressure
    nos = length(y); # number of substances
    P = Pestimate(temperature, Tc, Pc, Ac);
    P = 1 / sum(y / P);
    x = y / Kestimate(temperature, P, Tc, Pc, Ac);
    c = 0; # iteration control
    # 1. calc. a, b liquid phase
    abg = ABSRK(temperature,
                y,
                nos);
    repeat{
      vg = VOLSRK(phase=-1,    # 2. calc. v, phi liquid @ P
                  temperature,
                  P,
                  abg$a,
                  abg$b);
      phig = PHISRK(temperature,
                    P,
                    vg,
                    y,
                    abg$a,
                    abg$b,
                    abg$bi,
                    abg$aij,
                    nos=nos);
      abl = ABSRK(temperature, # 3. calc. a, b gas phase
                  x,
                  nos);
      vl = VOLSRK(phase=1,  # 4. calc. v, phi gas @ P
                  temperature,
                  P,
                  abl$a,
                  abl$b);
      phil = PHISRK(temperature, 
                    P,
                    vl,
                    x,
                    abl$a,
                    abl$b,
                    abl$bi,
                    abl$aij,
                    nos=nos); 
      x = y * phig / phil;  # 5. calc. K = phi_l / phi_g, S = sum(x * K)
      S = sum(x);
      if ((abs(S - 1) < e) || (c > 99)){break;}
      else{
        c = c + 1;
        S = 1;
        for(i in 1:nos){
          S = S - x[i];
        }
        temp = S * P;
        P = P + temp;
      }
    }
    result = list(pressure=P*1e+5,
                  x=x,
                  phil=phil,
                  phig=phig,
                  S=S,
                  iterations=c);
    return(result);
  } # end calc.dewpoint.P
  calc.dewpoint.T = function(P, y, ...){ # calculate dewpoint temperature
    P = P * 1e-5;
    nos = length(y); # number of substances
    temperature = Testimate(P, Tc, Pc, Ac);
    temperature = sum(y * temperature);
    x = y / Kestimate(temperature, P, Tc, Pc, Ac);
    x = x / sum(x)
    c = 0;
    repeat{
      abl = ABSRK(temperature, 
                  x,
                  nos);# 1. calc. a, b liquid phase
      vl = VOLSRK(phase=1,    # 2. calc. v, phi liquid @ P
                  temperature,
                  P,
                  abl$a,
                  abl$b);
      phil = PHISRK(temperature,
                    P,
                    vl,
                    x,
                    abl$a,
                    abl$b,
                    abl$bi,
                    abl$aij,
                    nos=nos);
      abg = ABSRK(temperature, # 3. calc. a, b gas phase
                  y,
                  nos);
      vg = VOLSRK(phase=-1,  # 4. calc. v, phi gas @ P
                  temperature,
                  P,
                  abg$a,
                  abg$b);
      phig = PHISRK(temperature, 
                    P,
                    vg,
                    y,
                    abg$a,
                    abg$b,
                    abg$bi,
                    abg$aij,
                    nos); 
      K = phil / phig;  # 5. calc. K = phi_l / phi_g, S = sum(x * K)
      x = y / K;
      S = sum(x);
      if ((abs(S - 1) < e) || (c > 99)){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (S - 1) / S;
        temperature = temperature + temp;
        x = x / S;
        x = x / sum(x);
      }
    }
    result = list(temperature=temperature,
                  x=x,
                  phil=phil,
                  phig=phig,
                  S=S,
                  iterations=c);
    return(result);
  } # end calc.dewpoint.T
  calc.saturation.P = function(temperature, x, ...){
    P = Pestimate(temperature, Tc, Pc, Ac);
    nos = 1; # number of substances
    c = 0; # iteration control
    abl = ABSRK(temperature, # 1. calc. a, b liquid phase
                x,
                nos);
    repeat{
      vl = VOLSRK(phase=1,    # 2. calc. v, phi liquid @ P
                  temperature,
                  P,
                  abl$a,
                  abl$b);
      phil = PHISRK(temperature,
                    P,
                    vl,
                    x,
                    abl$a,
                    abl$b,
                    abl$bi,
                    abl$aij,
                    nos=nos);
      abg = ABSRK(temperature, # 3. calc. a, b gas phase
                  y,
                  nos);
      vg = VOLSRK(phase=-1,  # 4. calc. v, phi gas @ P
                  temperature,
                  P,
                  abg$a,
                  abg$b);
      phig = PHISRK(temperature, 
                    P,
                    vg,
                    y,
                    abg$a,
                    abg$b,
                    abg$bi,
                    abg$aij,
                    nos=nos); 
      S = phil - phig;  # 5. calc. K = phi_l / phi_g, S = sum(x * K)
      if ((abs(S) < e) || (c > 10)){break;}
      else{
        c = c + 1;
        temp = P * S;
        P = P + temp;
      }
    }
    result = list(pressure=P*1e+5,
                  y=y,
                  phil=phil,
                  phig=phig,
                  S=S,
                  iterations=c);
    return(result);
  } # end calc.saturation.P
  calc.saturation.T = function(P, x, ...){
    P = P * 1e-5;
    temperature = Testimate(P, Tc, Pc, Ac);
    nos = length(x); # number of substances
    c = 0;
    repeat{
      abl = ABSRK(temperature, 
                  x=x,
                  nos=nos);# 1. calc. a, b liquid phase
      vl = VOLSRK(phase=1,    # 2. calc. v, phi liquid @ P
                  temperature,
                  P,
                  abl$a,
                  abl$b);
      phil = PHISRK(temperature,
                    P,
                    vl,
                    x,
                    abl$a,
                    abl$b,
                    abl$bi,
                    abl$aij,
                    nos=nos);
      abg = ABSRK(temperature, # 3. calc. a, b gas phase
                  y,
                  nos);
      vg = VOLSRK(phase=-1,  # 4. calc. v, phi gas @ P
                  temperature,
                  P,
                  abg$a,
                  abg$b);
      phig = PHISRK(temperature, 
                    P,
                    vg,
                    y,
                    abg$a,
                    abg$b,
                    abg$bi,
                    abg$aij,
                    nos=nos); 
      S = phig / phil
      if ((abs(S - 1) < 1e-4) || (c > 199)){break;}
      else{
        c = c + 1;
        temp = 0.1 * temperature * (S - 1) / S;
        temperature = temperature + temp;
      }
    }
    result = list(temperature=temperature,
                  y=y,
                  phil=phil,
                  phig=phig,
                  S=S,
                  iterations=c);
    return(result);
  } # end calc.saturation.T
  if((is.null(P)) && (!is.null(temperature)) && (!is.null(x)) && (length(x) > 1)){ # bubblepoint pressure = f(T,x)
    return(calc.bubblepoint.P(temperature, x));
  }
  if((is.null(P)) && (!is.null(temperature)) && (is.null(x)) && (!is.null(y))){ # dewpoint pressure = f(T,y)
    return(calc.dewpoint.P(temperature, y))  
  }
  if((is.null(P)) && (!is.null(temperature)) && (is.null(x)) && (is.null(y))){
     warning('No x or y given to calculate bubble- (f(x)) or dewpoint pressure (f(y)) at ', temperature, 'K.');
  }
  if((is.null(temperature)) && (!is.null(P)) && (!is.null(x)) && (length(x) > 1)){ # bubblepoint temperature = f(P,x)
    return(calc.bubblepoint.T(P, x));
  }
  if((is.null(temperature)) && (!is.null(P)) && (is.null(x)) && (!is.null(y))){ # dewpoint temperature = f(P,y)
    return(calc.dewpoint.T(P, y));
  }
  if((is.null(temperature)) && (!is.null(P)) && (is.null(x)) && (is.null(y))){
    warning('No x or y given to calculate bubble- (f(x)) or dewpoint temperature (f(y)) at ', P, ' bar.');
  }
  if((length(x) == 1) && (!is.null(temperature)) && (is.null(P)) && (is.null(y))){ # saturation pressure calculation
    y = 1;
    return(calc.saturation.P(temperature, x));    
  }
  if((length(x) == 1) && (is.null(temperature)) && (!is.null(P)) && (is.null(y))){ # saturation temperature calculation
    y = 1;
    return(calc.saturation.T(P, x));
  }
}