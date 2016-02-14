Substance = function(name){
  AntoineVals = read.csv('Substances.Antoine-tbl.csv', sep=','); # load Antoine coefficents table
  Antoine = AntoineVals[grep(paste('^', name, sep='', ' '), AntoineVals$Substance),];
  if(nrow(Antoine) == 0){
    Antoine = AntoineVals[grep(paste('^', name, sep='', '$'), AntoineVals$Substance),];
  }
  if(nrow(Antoine) == 0){
    stop('\nAntoine coefficients for ', name,' not found.');
  }
  Antoine$Substance = NULL;
  criticalVals = read.csv('Substances.critical-tbl.csv', sep=','); # load critical data table
  criticals = criticalVals[grep(paste('^', name, sep='', ' '), criticalVals$Substance), ];
  if(nrow(criticals) == 0){
    criticals = criticalVals[grep(name, criticalVals$Substance), ];
  }
  if(nrow(criticals) == 0){
    stop('\nCritical values for ', name,' not found.');
  }
  EleMass = read.csv('ElementMasses.csv', sep=',');  # load element masses file to calc. the molare mass
  tmp1 = unlist(strsplit(gsub('[[:digit:]]', ' ',criticals$Formula), ' ')) # get characters
  tmp1 = tmp1[tmp1 != '']; # remove empty strings
  tmp2 = as.numeric(unlist(strsplit(gsub('[[:alpha:]]', ' ',criticals$Formula), ' '))); # get numbers
  tmp2 = tmp2[!is.na(tmp2)]; # remove NA
  M = 0; # molare mass
  for (i in 1:length(tmp1)){
    mass = EleMass[grep(tmp1[i], EleMass$Element), ]$Mass;
    if(i <= length(tmp2)){
      M = M + tmp2[i] * mass;
    }
    else{
      M = M + mass;
    }
  }
  string = tail(unlist(strsplit(gsub('',' ', criticals$Formula), ' ')), -1);
  elements = string[1];
  mass = NULL;
  factors = NULL;
  # 2. start from index two (first is upper letter) and check for upper letter
  i = 2
  while (i <= length(string)){
    if(grepl('[A-Z]', string[i])){ # is upper?
      elements = append(elements, string[i]);
      if(!grepl('[0-9]', string[i - 1])){  
        factors = append(factors, 1);
      }
      if(i == length(string)){
        factors = append(factors, 1);
      }
    }
    else if(grepl('[a-z]', string[i])){ # is lower?
      elements[i - 1] = paste(elements[i - 1], string[i], sep='');
    }
    else{ # is number!
      tmp = '';
      while(grepl('[0-9]', string[i])){
        tmp = paste(tmp, string[i], sep='')
        i = i + 1;
      }
      i = i - 1;
      factors = append(factors, as.numeric(tmp));
    }
    i = i + 1;
  }
  for (i in 1:length(elements)){
    mass = append(mass, EleMass[grep(elements[i], EleMass$Element), ]$Mass);
  }
  M2 = sum(factors * mass);
  result = list(Antoine = Antoine,
                Ac = criticals$Ac,
                Pc = criticals$Pc,
                Tc = criticals$Tc,
                MolarMass1 = M,
                MolarMass2 = M2);
  return(result)
}
# Substances=NULL;
# for(i in 1:2){
#   Substances=list(Substances, Substance(subs[i]));
# }
# names(Substances) = name;