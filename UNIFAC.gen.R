UNIFAC.gen = function(Substances=NULL,
                      verbose=FALSE){
  # tool to sum the UNIFAC parameters and to generate the surface and
  # the interaction parameters matrizes
  # 
  # input: vector with substances or NULL
  #
  # return: list with 2 matrices: $unu = surface parameters
  #                               $aij = interaction parameters
  # file: UNIFAC-tbl3.csv contains example substances with group contributions
  #       all new substances are stored in it
  # define functions:
  # 'printf' Formated printing, example: printf('Something equals %d',1)
  printf = function(...) cat(sprintf(...));
  # 'get.aij'
  get.aij = function(i,j){# read UNIFAC interaction paramters from csv file
    atbl = read.csv('UNIFAC-tbl2.csv',sep=',');
    # search insid the UNIFAC interaction parameter table
    for (k in 1:nrow(atbl)){
      if ((atbl[k,1] == i) && (atbl[k,2]== j)){
        return(atbl[k,3]);
        break;
      }
      else if (k == nrow(atbl)){
        printf('Combination of group %d and group %d not found.',i,j);
        return(0);
      }
    }
  }
  # 'get.UNIFAC' Get UNIFAC paramters from table get.unifac('formula') 
  # -> k, Group, Rk, Qk
  get.UNIFAC = function(formula){
    utbl = read.csv('UNIFAC-tbl1.csv',sep=','); # read UNIFAC surface paramters
    # search the UNIFAC table for functional group by formula and 
    # return the line or 0
    # formula = case sensitive string e.g. H2O, OH, CH2=CH, Morpholine
    line = utbl[grep(paste('^',formula,sep='',' '),utbl$Group),];
    if (nrow(line)==0){
      line = utbl[grep(paste('^',formula,sep='','$'),utbl$Group),];
    }
    if (nrow(line)==0){
      stop('\nGroup ',formula,' not found.\n');
    }
    return(line);
  }
  # 'get.sub' Read the functional group contribution, by name
  get.sub = function(name){
    subtbl = read.csv('UNIFAC-tbl3.csv',sep=',');
    sub = subtbl[grep(paste('^',name,sep='',' '),subtbl$Name),];
    if (nrow(sub)==0){
      sub = subtbl[grep(paste('^',name,sep='','$'),subtbl$Name),];
    }
    if (nrow(sub)==0){
      printf('\nSubstance "%s" not found.\n',name);
    }
    else{
      return(sub)
    }
  }
  # define variables
  sn = Substances;
  k = c(); # subgroup number
  if(verbose){
      printf('Welcome to the "UNIFAC MATRIZES GENERATOR"\n\n');
      printf('All known substances are stored in file "UNIFAC-tbl3.csv",\n');
      printf('new ones can be defined by name and group contribution.');
  }
  if(is.null(sn)){
    printf('Generate the UNIFAC Matrix in two steps:\n\n');
    printf('1.Type the names of all substances.\n')
    printf('2. Type the group contribution for unknown substances.\n\n');
    printf('Type the names of all substances (case sensitive)!');
    # 1. read the names from user input
    sn = readline(); # substances names
    sn = strsplit(sn,',')[[1]]; # split names at seperator ','
  }
  if(verbose){
    printf('The following substances are in use:\t');
    printf('%s\t', sn);
  }
  nos = length(sn); # number of substances
  tmp1 = data.frame();
  for (i in 1:nos){
    tmp2 = get.sub(sn[i]); # check for known substances in file 'UNIFAC-tbl3.csv'
    if(!is.null(tmp2)){
      tmp1 = rbind(tmp1,tmp2); # one substance per row
    }
    else{ # define substance
      tmp2 = data.frame(Name=sn[i], Groups=NA)[, ]
      printf('\nEnter the group contribution for %s:\t[1*CH3 1*CH2 1*OH]',sn[i]);
      tmp2[2] = readline();
      tmp1 = rbind(tmp1,tmp2);
      write.table(tmp2,file='UNIFAC-tbl3.csv',sep=',',
                  row.names=F,col.names=F,append=T);
    }
  } # transform listed functional groups in matrix
  groups = strsplit(paste(tmp1[,2]),'\\t'); # 1st split -> groups per substances
  tmp1 = data.frame();
  for (i in 1:length(groups)){ # from 1st to last substance
    groups[i] = strsplit(groups[[i]],'\\ '); # 2nd split between groups
    nu = rep(0,1);
    gr = rep(NA,1);
    for (j in 1:length(groups[[i]])){ # from 1st to last group per substance
      tmp1 = strsplit(groups[[i]][j],'\\*')[[1]]; # 3rd split -> amount of groups
      nu[j] = tmp1[1];
      gr[j] = tmp1[2];
    }
    if (i == 1){ # 1st substance
      num = cbind(gr,nu); # matrix of group and amount per group
    }
    else{
      tmp2 = rep(0,nrow(num)); # zero column
      num = cbind(num,tmp2); # new substance -> new column
      for(l in 1:length(gr)){ # from 1st to last group per substance
        pos = grep(paste('^',gr[l],sep='','$'),num[,1]); # group exists?
        if (length(pos) > 0){
          num[pos,i+1] = nu[l]; # hit -> set value
        }
        else{ # new group
          tmp2 = rep(0,ncol(num)-2); # repeat zeros for the new line
          tmp3 = c(gr[l],tmp2,nu[l]); # vector with group, 0, amount of group
          num = rbind(num,tmp3); # add the new line
        }
      }
    }
  }
  nog = nrow(num); # number of groups
  tmp1 = data.frame();
  for (i in 1:nog){
    tmp2 = get.UNIFAC(num[i,1]); # get the UNIFAC parameters
    tmp1 = rbind(tmp1,tmp2);
  }
  if((nrow(tmp1) == nog) || (length(tmp1) == nog) ){
    unu = cbind(tmp1,num[,-1]);
  }# data frame with UNIFAC values and amount per substance
  else{
    stop('Not all groups have been identified. 
         May check for group with get.UNIFAC(group)')
  }
  for(i in 5:ncol(unu)) unu[,i] = as.numeric(as.character(unu[,i]));
  colnames(unu)[5:ncol(unu)] = sn;
#  unu = data.frame(unu, row.names = NULL);
  # 3. Generate the interaction parameters matrix for the subgroups
  aij = matrix(0,nog,nog);
  for (i in 1:nog){
    for (j in 1:nog){
      aij[j,i] = get.aij(unu[j,2],unu[i,2]);
    }
  }
  if(verbose){
    printf('\n\nThe UNIFAC matrix is ready for use.\n\n');
    print(unu);
    printf('\n\nThe interaction parameters matrix is ready for use.\n\n');
    print(aij);
  }
  result = list(unu=unu,
                aij=aij);
  return(result);
}