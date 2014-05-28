# load functions
source('compcheck.R')
source('Antoine.R')
source('Tboil.R')
source('unifac.R')
Txvle2 <- function (compone,comptwo,unu,aij,P) # main function, co functions are below
{
  #  Author:  Hannes Seidel  Hannes.Seidel@tu-cottbus.de   
  #        			                                          
  #  Name:	Txvle		Version:	0.2 			                  
  #
  #  This function calculate the molar fractions of the     
  #  liquid and gas phase for binary mixtures.              
  #  It starts at the boiling point of the highboiler.       
  #  The saturated pressure and the activity coefficient    
  #  are calculated using the Antoine equation and the      
  #  UNIFAC model.                                     
  #  The gas phase is ideal.                                
  #                                                         
  #  Txvle (lowboiler,highboiler,u,r,q,system pressure)     
  #                                                         
  #  u = [2,2] matrix                                       
  #  r,q = 2 dim vector                                     
  #                                                         
  #  returns a resultant matrix with T,x,y,act1,act2,setps
  
  ## x -> psat, act @ T.boiling -> sum (x * psat * act) == P ? -> T(P,x), y(P,T,x)
  
  comp1<-compcheck(compone) # load Antoine constants
  comp2<-compcheck(comptwo)
  e <- 0.001 # best assumption, 0.0001 creates failures
  T0 <- Tboil(comptwo,P) # estimate lower boiling point
  x<-seq(0,1,0.01) # x vector
  txy<-c() # initiate the resultant vector
  Temperature <- T0 # start point (ideal)
  for (i in 1:length(x))
  {
    c <- 0 # iteration control (emergency break)
    repeat
    {
      Ps1 <- Antoine(comp1[2],comp1[3],comp1[4],Temperature) # saturated pressure
      Ps2 <- Antoine(comp2[2],comp2[3],comp2[4],Temperature)
      act <- unifac(c(x[i],1-x[i]),unu,aij,Temperature) # activity coefficients
      Pc <- (x[i]*Ps1*act[1] + (1-x[i])*Ps2*act[2])
      y <- c(x[i]*act[1]*Ps1/Pc, (1-x[i])*act[2]*Ps2/Pc) # gas fraction y1 <- x1*Ps1*act1 / P
      rat <- Pc / P # Pcalc / P - ratio
      if ((abs(rat - 1) < e)|| (c > 999)) {break} # solution criteria and emergency break
      else if (rat < 1) {Temperature <- Temperature + rat} # decrease temperature with the rat
      else {Temperature <- Temperature - rat} # increase
      c <- c + 1 # count the iteration per x step
    }
    txy <- c(txy,Temperature,x[i],y[1],act,c,Pc) # fill the resultant vector
  }
  TXY<-matrix(txy,nrow=length(x),ncol=7,byrow=T) # transform to resultant matrix
  colnames(TXY)<-c("T [K]","x(T,P)","y(T,P,x)","act 1","act 2","steps","Pc") # set the column names
  return(TXY)
}

## Examples, uncomment for using
# Example EtOH-H2O, # EtOH: compvals[5,], H20: compvals[30,], P = Patm = 101325 Pa
# ulist <- unifactool() # interactive user input
# unu <- ulist[[1]] # UNIFAC surface parameters and component matrix
# aij <- ulist[[2]] # UNIFAC interaction parameters matrix [k,k]
#TXY2 <- Txvle2(5,30,unu,aij,1.013e+5)
