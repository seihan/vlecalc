solvemixtureproblem = function(substances = NULL,
                               pressure = NULL,
                               boilline = NULL){
  source('os.boilingline.R')
  objF = function(fractions, ...){
    fractions = abs(fractions);
    fractions = fractions / sum(fractions);
    blp = os.boilingline(substances, fractions, pressure);
    diff = sum(abs(polyval(p=blp, boilline[,1]) - boilline[,2]));
    print(diff)
    x = seq(0, 100, 1);
    lines(x, polyval(blp, x), col='red');
    return(diff);
  }
  plot(boilline[,1], boilline[,2]);
  nos = length(substances);
  fractions = rep(1, nos) / nos;
  result = optim(par=fractions, fn=objF);
  return(result);
}