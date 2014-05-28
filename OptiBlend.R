OptiBlend <- function(){
  # Author:  Hannes Seidel  Hannes.Seidel@tu-cottbus.de   
  #            	                                          
  # Name:	optiblend		Version:	0.1 			                  
  #
  # Optimize surrogate mixtures to reach vaporization lines of given fuels
  #
  # 1. set initial mixture: x_i = 1 / number of components
  # 2. simulate the vaporization line and compare it to a measured one
  # 3. reblend mixture by minimizing the difference
  #
  # f = dT(bubble) + dT(dew)
  # min{f}
  
}
test6<-Siedelinie.mc.real(c(30,24,23,14,17,15,27,13,16),c(1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9),unu,aij,1.013e+5)
plot(100-test6[[5]][,2],test6[[5]][,1],type='l',lwd=2,main='Vaporization Line - FACE#6',xlab='Evaporated fraction [%]',ylab='T [K]')
