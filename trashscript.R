# temperature <- readline('Enter initial guess for T (K).\n[Default is to use last value.\nProgram will make an initial guess if you enter zero or leave it blank]\n ')
# temperature <- as.numeric(temperature)
# 
# if(is.na(temperature)||(temperature<=0))
# {
#   print('No Temperature, no problem.')
#   print('9')
# }else
# {
#   print(temperature)
# }
# pdf("~/Dropbox/MA-Thesis/Fuels/calculations/xy-etoh-h2o_laar-margules-srk-uniquac-ddb.pdf")
# plot(c(0,1),c(0,1),type="l",main="Vapor-Liquid-Equilibrium Ethanol / Water \nat 300 K",xlab=expression("x"[Ethanol]*" [mol/mol]"),ylab=expression("y"[Ethanol]*" [mol/mol]"))
# lines(exvals[,1],exvals[,2],type='p',pch=8,lwd=1)
# lines(ilmXYP[,1],ilmXYP[,2],type='l',lty=1,lwd=2,col='gray66')
# lines(ilmXYP[,1],ilmXYP[,3],type='l',lty=6,lwd=2.5,col='blue2')
# lines(ilmXYP[,1],ilmXYP[,4],type='l',lty=4,lwd=2.5,col='red')
# lines(srkXYP[,1],srkXYP[,2],type='l',lty=2,lwd=2.5,col='firebrick')
# lines(uXYP[,1],uXYP[,2],type='l',lty=2,lwd=2.5,col='seagreen')
# legend(0.05,1,c("Experiment (DDB)","Ideal","Van Laar","Margules","SRK (k12 = -0.1)","UNIQUAC (PHI = 1)"),lty=c(NA,1,6,4,2,2),lwd=c(1.5,2,2,2,2),pch=c(8,NA,NA,NA,NA,NA),col=c("black","gray66","blue2","red","firebrick","seagreen"),cex=0.7,merge=FALSE)
# dev.off()

# ln phi0j (T,P) = Bjj(T) * P / R*T
# Tsonopoulos
# Bjj * Pc / R * Tc = f0 + f1*Wc
# f0 = 1.445 - 0.33 / Tr - 0.1385 / Tr² - 0.0121 / Tr³ - 0.000607 / Tr^8
# f1 = 0.0637 + 0.331 / Tr² -0.423 / Tr³ - 0.008 / Tr^8

Tsonopoulos <- function(Pc,Tc,Wc,Temperature)
{
  R <- 83.1433
  Tr <- Temperature / Tc
  f0 = 0.1445 - (0.33 / Tr) - (0.1385 / Tr^2) - (0.0121 / Tr^3) - (0.000607 / Tr^8)
  f1 = 0.0637 + (0.331 / Tr^2) - (0.423 / Tr^3) - (0.008 / Tr^8)
  B <- ((R*Tc)/Pc)*(f0 + f1*Wc)
  return(B)
}
Bmix1 <- function(B11,B22)
{
  B12 <- sqrt(B11*B22)
}
Bmix2 <- function(B11,B22)
{
  B12 <- (1/8)*((B11)^(1/3)+(B22)^(1/3))^3
}
viriphi <- function(B,P,Temperature)
{
  R <- 83.1433
  lnphi <- B*P/(R*Temperature)
  return(exp(lnphi))
}
#xydioxwater<-read.csv('~/Dropbox/MA-Thesis/Fuels/calculations/values/exvals_dioxan-wasser.csv',head=T,sep=',')
EtOH <- lookupcritical(17)
Tc1 <- as.numeric(EtOH[2])
Pc1 <- as.numeric(EtOH[3])
Wc1 <- as.numeric(EtOH[4])
H2O <- lookupcritical(16)
Tc2 <- as.numeric(H2O[2])
Pc2 <- as.numeric(H2O[3])
Wc2 <- as.numeric(H2O[4])
comp1 <- compcheck(4)
comp2 <- compcheck(9)
Temperature <- 354.0114
x <- 0.5
e <- 0.001
c <- 0
repeat
{
  Ps1 <- Antoine(comp1[2],comp1[3],comp1[4],Temperature)
  Ps2 <- Antoine(comp2[2],comp2[3],comp2[4],Temperature)
  B1 <- Tsonopoulos(Pc1,Tc1,Wc1,354.0114)#Temperature)
  B2 <- Tsonopoulos(Pc2,Tc2,Wc2,354.0114)#Temperature)
  phi1 <- viriphi(B1,patm,Temperature)
  phi2 <- viriphi(B2,patm,Temperature)
  rat <- P / (Ps2/(1-x*(Ps1 - Ps2)/Ps1))
  if((abs(rat - 1) < e)||(j > 10000)){break}
  else if (rat < 1){T100 <- T100 - rat}
  else {T100 <- T100 + rat}
  j <- j + 1
  if ((y < e) || (c > 999)){break}
  else {Temperature <- Temperature - 0.1}
  c <- c + 1
}
## gnomesort
# Sort vector
gnomesort <- function(vec){
  i <- 1  
  while (i <= length(vec)){
    if (i == 1 || vec[i-1] <= vec[i]){
      i <- i + 1
    }
    else{
      tmp <- vec[i]
      vec[i] <- vec[i-1]
      i <- i - 1
      vec[i] <- tmp
    }
  }
  return(vec)
}
## gnomesort.matrix
# Sort matrix by specific column
gnomesort.matrix <- function(mat,col){
  i <- 1  
  while (i <= nrow(mat)){
    if (i == 1 || mat[i-1,col] <= mat[i,col]){
      i <- i + 1
    }
    else{
      tmp <- mat[i,]
      mat[i,] <- mat[i-1,]
      i <- i - 1
      mat[i,] <- tmp
    }
  }
  return(mat)
}

# multicore
library(foreach)
library(doMC)
library(multicore)
ncore = multicore:::detectCores()
registerDoMC(cores = ncore)
results <- foreach(i = 1:5, .combine=c) %dopar% {
  i+i
}
results
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
