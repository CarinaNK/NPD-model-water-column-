#marine modelling NPD

rm(list=ls())
library(deSolve)
library(spam)
library(viridis)
library(viridisLite)
library(fields)


NPD = function(dz,Hn){
  param=c()
  param$dz=dz  
  param$z=seq(param$dz/2,200,by=param$dz)
  n <- length(param$z)
  param$kp=0.1 #m^2/my mol N, self shading coeff.
  param$kw=0.1#m^-1 light absorbed by water.
  param$I0=350  # mymol photons*m^2/sec
  param$Hi=30 # maybe change the unit? mymol photons/m2/s halfsaturation of light
  param$Hn=Hn# mymol/ m^3
  param$m=0.24#mortality, pr day
  gmax=1.5 # pr dag
  param$u=3# this is dropping of plankton parm. m/day
  param$D=5*10^-5*3600*24 #m^2/day diffusivity.
  param$nbot=30
  param$alpha=1
  param$eps=0.1 # day ^-1, remineralization rate.
  param$bet=5 # strength of light effect, just a variable skalar.
  param$gamma=0.2 # grazing m^3 (mmol N)^-1 day^-1
  
  y=c(1:(n*3))*0
  y[1]=1 # initial phytoplankton.
  y[(n+1):(2*n)]=20# initial nutrients
  y[(2*n+1):(3*n)]=5# initial detritus.
  
  lightemissions=function(P,Detrit,param,t){
    summerpleb=(cos(2*pi*(t/365))+1.01)*param$bet
    Q=param$kp*param$dz*((cumsum(P)-P/2)+(cumsum(Detrit)-Detrit/2))
    I=param$I0*exp(-param$kw*param$z-Q)*summerpleb
    return(I)
  }
  
  
  Ja=c() # advective flux
  Jd=c() #diffusive flux
  Jdn=c(0) # nutrient flux
  Jdt=c()  #detritus flux
  Jadt=c() # also detritus flux
  
  P=c()
  N=c()
  dPdt=c()
  dNdt=c()
  
  
  diffusivefux=function(t,y,mu){
    
    P=y[1:n]
    N=y[(n+1):(n*2)]
    Detrit=y[(2*n+1):(3*n)]
    for (i in 2:n){
      Ja[i]=param$u*P[i-1]
      Jd[i]=-param$D*((P[i]-P[i-1])/param$dz)
      Jdn[i]=-param$D*((N[i]-N[i-1])/param$dz)
      Jadt[i]=param$u*Detrit[i-1]
      Jdt[i]=-param$D*(Detrit[i]-Detrit[i-1])/param$dz
      
    }
    Ja[c(1,n+1)]=0  
    Jd[c(1,n+1)]=0
    Jdn[c(1,n+1)]=c(0, -param$D*(param$nbot-N[n])/param$dz) 
    Jadt[c(1,n+1)]=c(0,param$u*Detrit[n])
    Jdt[c(1,n+1)]=0  
    J=Ja+Jd
    JD=Jadt+Jdt
    
    
    ii=lightemissions(P,Detrit, param,t)

    g=gmax*pmin(param$Hi*ii /sqrt(gmax^2 + (param$Hi*ii)^2),N/(N+param$Hn))
    
    dPdt=-((J[2:(n+1)]-J[1:(n)])/param$dz) + (g*P)-param$m*P-param$gamma*P*P
    dNdt=-g*P+param$eps*Detrit-((Jdn[2:(n+1)]-Jdn[1:(n)])/param$dz) 
    dDdt=-((JD[2:(n+1)]-JD[1:(n)])/param$dz)+param$m*P+param$gamma*P*P-param$eps*Detrit
    
    
    Y=c(dPdt,dNdt,dDdt)
    list(Y)
  }
  
  timey=3000
  
  res=ode.1D(y=y,func=diffusivefux, times=0:timey,parms=1,nspec=1)
  
  return(res)


}

res =NPD(1,0.3)
restid = res[,1]

z=seq(1/2,200,by=1)
n <- length(z)
resP=res[,2:(n+1)]
resN=res[,(n+2):(2*n+1)]
resD=res[,(2*n+2):(3*n+1)]

image.plot(x=restid,y=z,z=resP,main="Phytoplankton [mmol/m^3]",ylim = rev(range(z)),xlab="Time [d]",ylab="Depth [m]")
image.plot(x=restid,y=z,z=resN,main="Nutrients [mmol/m^3]",ylim = rev(range(z)),xlab="Time [d]",ylab="Depth [m]")
image.plot(x=restid,y=z,z=resD,main="Detritus [mmol/m^3]",ylim = rev(range(z)),xlab="Time [d]",ylab="Depth [m]")



# sensitivity analysis dz
res = NPD(0.2,0.3)
z=seq(0.2/2,200,by=0.2)
n <- length(z)
resP=res[,2:(n+1)]
timey=3000


plot(resP[(timey+1),], z, ylim = rev(range(z)), type='l', lwd=0.25, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]")

dz = c(0.4,0.8,1.6,3.2,6.4)
lwd = c(0.5,1,2,3,4)
c=1

for (i in dz){
  res = NPD(i,0.3)
  z=seq(i/2,200,by=i)
  n <- length(z)
  resP=res[,2:(n+1)]
  lines(resP[(timey+1),], z, ylim = rev(range(z)), pch = 19, lwd =lwd[c])
  c=c+1
}



# sensitivity analysis Hn
res = NPD(1,0.3)
z=seq(1/2,200,by=1)
n <- length(z)
resP=res[,2:(n+1)]
timey=3000


plot(resP[(timey+1),], z, ylim = rev(range(z)), type='l', lwd=1, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]")

Hn = c(0.01,0.1,1,10)
lwd = c(0.25,0.5,2,3)
c=1

for (i in Hn){
  res = NPD(1,i)
  resP=res[,2:(n+1)]
  lines(resP[(timey+1),], z, ylim = rev(range(z)), pch = 19, lwd =lwd[c])
  c=c+1
}

legend(0.35, 0, legend=c("Hn = 0.01", "Hn = 0.1", "Hn = 0.3", "Hn = 1", "Hn = 10"),lwd= c(0.25,0.5,1,2,3), cex=0.8)

