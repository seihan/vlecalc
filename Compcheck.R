Compcheck<-function(comp) {
  # returns a line with Antoine constants
  # comp must be a number (for now, improvements will follow)
  compvals <- read.csv('comvals.csv',sep=',')
  return(compvals[comp, ])
}