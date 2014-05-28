unifactool <- function(){
  source('printf.R')
  source('get.unifac.R')
  source('get.aij.R')
  # tool for sum the unifac parameters
  # it returns a list where list [[1]] contains the surface-, [[2]] contains the interaction parameters
  noc <- 1 # number of components
  k <- c() # subgroup number
  printf('Welcome.\nGenerate the UNIFAC Matrix in two steps:\n\n1.\tType functional groups of all substances')
  printf(' in\n\t\t\tthe right order.\n2. Tell the amount of each group in substance.\n\n')
  printf('1st step.\nType the groups of all substances (case sensitive):\t[CH3,CH2,SiH2,...]')
  # 1. read the names
  sgn <- readline() # subgroup name
  sgn<-strsplit(sgn,',')[[1]]
  nog <- length(sgn) # number of groups
  # 2. read the amount of group per substance
  nu <- rep(0,nog) # vector with amount of groups per substance (nu)
  printf('\n2nd step.\n')
  repeat{
    for(i in 1:nog){
      printf('Amount of group %s in substance %d ?',sgn[i],noc)
      nu[i] <- as.numeric(readline())
    }
    if(noc == 1){
      num<-data.frame('nu 1'= nu) # matrix with number of subgroups for all components
    }
    else{
      num[,paste(noc)] <- nu
    }
    sv <- readline('You like to add another substance?\t[Y/n]\t\n')
    tmp1 <-data.frame()
    for (i in 1:nog){
      tmp2 <- get.unifac(sgn[i])
      tmp1 <- rbind(tmp1,tmp2)
    }
    if((nrow(tmp1)==nrow(num))||length(tmp1)==length(num)){
      unu <- cbind(tmp1,num) # create data frame with UNIFAC values and amount per substance
    }
    else{
      stop('Something wrong. May check for group with get.unifac(group)')
    }
    if ((sv == 0) || (sv == 'n')){break}
    else{
      noc <- noc + 1 # score the number of components
    }
    # 3. Generate the interaction parameters matrix for the subgroups
    aij <- matrix(0,nog,nog)
    for (i in 1:nog){
      for (j in 1:nog){
        aij[j,i] <- get.aij(unu[j,1],unu[i,1])
      }
    }
  }
  return(list(unu,aij))
}