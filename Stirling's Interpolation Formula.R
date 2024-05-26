# Stirling's Interpolation Formula
Stirling<-function(x, y, x0){
  h<-x[2]-x[1]
  n<-length(which(x<x0))-1
  u<-(x0-x[n+1])/h
  v<-rep(0,n)
  v[1]<-u^2
  for(i in 2:n){
    v[i]<-(u^2-i+1)*v[i-1]
  }
  S<-0
  for(i in 1:n){
    del1<-diff.default(y, lag=1, differences=2*i-1)
    del2<-diff.default(y, lag=1, differences=2*i)
    S<-S+((v[i]/(u*factorial(2*i-1)))*(del1[n+1-i]+del1[n+2-i])/2    +(v[i]/(factorial(2*i)))*(del2[n+1-i]))
  }
  return(y[n+1]+S)
}
# An example
x<-seq(.51, .57, .01)
y<-c(.5292437, .5378987, .5464641, .5549392, .5633233, .5716157, .5798158)
Stirling(x, y, 0.5437)	# More appropriate than Bessel’s
# R output
[1] 	0.5580519606
Stirling(x, y, 0.55)	#Bessel’s is more appropriate
# R output
[1] 	0.5633233
