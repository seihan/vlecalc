Unifactool2 <- function(){
  source('printf.R')
  source('get.unifac.R')
  source('get.aij.R')
  source('get.sub.R')
  # tool for sum the unifac parameters
  # it returns a list where list [[1]] contains the surface-, [[2]] contains the interaction parameters
  k <- c() # subgroup number
  printf('Welcome.\nGenerate the UNIFAC Matrix in two steps:\n\n1.\tType the names of all substances.')
  printf('\n2. Type the group contribution for unknown substances.\n\n')
  printf('1st step.\nType the names of all substances (case sensitive):\t[n-decane,ethanol,water,...]')
  # 1. read the names from user input
  sn <- readline() # substances names
  sn<-strsplit(sn,',')[[1]] # split names at seperator ','
  nos <- length(sn) # number of substances
  tmp1 <- data.frame()
  for (i in 1:nos){
    if(nrow(get.sub(sn[i]))!=0){
      tmp2 <- get.sub(sn[i]) # check for known substances in file 'unifac-tbl3.csv'
      tmp1 <- rbind(tmp1,tmp2) # one substance per row
    }
    else{ # define substance
      tmp2 <- data.frame(Name=sn[i], Groups=NA)[, ] # names in frames have to be equal 
      printf('\nEnter the group contribution for %s:\t[1*CH3 1*CH2 1*OH]',sn[i])
      tmp2[2]<- readline()
      tmp1 <- rbind(tmp1,tmp2)
      write.table(tmp2,file='unifac-tbl3.csv',sep=',',row.names=F,col.names=F,append=T)
    }
  } # transform listed functional groups in matrix
  groups <- strsplit(paste(tmp1[,2]),'\\t') # 1st split between groups per substances
  tmp1 <- data.frame()
  for (i in 1:length(groups)){ # from 1st to last substance
    groups[i] <- strsplit(groups[[i]],'\\ ') # 2nd split between groups
    nu <- rep(0,1)
    gr <- rep(NA,1)
    for (j in 1:length(groups[[i]])){ # from 1st to last group per substance
      tmp1 <- strsplit(groups[[i]][j],'\\*')[[1]] # 3rd split between amount of groups
      nu[j] <- tmp1[1]
      gr[j] <- tmp1[2]
    }
    if (i == 1){ # 1st substance
      num <- cbind(gr,nu) # matrix of group and amount per group
    }
    else{
      tmp2 <- rep(0,nrow(num)) # zero column
      num <- cbind(num,tmp2) # new substance -> new column
      for(l in 1:length(gr)){ # from 1st to last group per substance
        pos <- grep(paste('^',gr[l],sep='','$'),num[,1]) # check for existing group
        if (length(pos) > 0){
          num[pos,i+1] <- nu[l] # hit -> set value
        }
        else{ # new group
          tmp2 <- rep(0,ncol(num)-2) # repeat zeros for the new line
          tmp3 <- c(gr[l],tmp2,nu[l]) # vector with group, zeros, amount of group
          num <- rbind(num,tmp3) # add the new line
        }
      }
    }
  }
  nog <- nrow(num) # number of groups
  tmp1 <-data.frame()
  for (i in 1:nrow(num)){
    tmp2 <- get.unifac(num[i,1]) # get the UNIFAC parameters from file 'unifac-tbl1.csv'
    tmp1 <- rbind(tmp1,tmp2)
  }
  if((nrow(tmp1) == nog) || (length(tmp1) == nog) ){
     unu <- cbind(tmp1,num[,-1]) # create data frame with UNIFAC values and amount per substance
  }
   else{
     stop('Not all groups have been identified. May check for group with get.unifac(group)')
   }
  
  names(unu)[5:length(unu)]<-sn # store the names of substances as colnames
  # 3. Generate the interaction parameters matrix for the subgroups
  aij <- matrix(0,nog,nog)
  for (i in 1:nog){
    for (j in 1:nog){
      aij[j,i] <- get.aij(unu[j,1],unu[i,1])
    }
  }
  return(list(data.matrix(unu),aij))
}