#INLA for Singapore/New Zealand data
source('read_X.r')

w=dlnorm(seq(1,50),1.132,0.742)
w=w/sum(w)

library(INLA)
start=322; end=nrow(X_sd) #year 2021, Singapore
start=552 #for NZ 
y=I[start:end]
Infect=rep(0,end-start+1)
for(i in start:end)
{
  if(i>50){
    Infect[i-start+1]=sum(rev(w)*I[(i-50):(i-1)])
  }else{
    Infect[i-start+1]=sum(rev(w[1:(i-1)])*I[1:(i-1)])
  }
}

lag=0
data=data.frame(cbind(y,Infect,X_sd[start:end-lag,]))
data$idx = 1:nrow(data)
f1 <- as.formula(paste0("y ~ 1 + f(idx, model='iid') + ", paste0(colnames(X_sd), collapse = " + ")))
ll=list()
for(i in 1:ncol(X_sd))
{
  ll[[colnames(X_sd)[i]]]=X_sd[start:end-lag,i]
}
ll[["(Intercept)"]]=rep(1,end-start+1)
lc=inla.make.lincombs(ll)
model=inla(f1,data=data,family = "poisson",control.family=list(link='log'), E=Infect, lincomb =lc,control.compute=list(dic=TRUE, cpo=TRUE, return.marginals=TRUE, return.marginals.predictor=TRUE, config = TRUE),control.predictor = list(compute=TRUE))
figure=FALSE
setwd(paste0(address,'data/',region,'/INLA'))
source(paste0(address,'coding/EpiMix_INLA_plot.r'))
model$dic$dic #DIC of the model

#obtain posterior \mu_t (expectation of R_t) estimates
Mu=R=matrix(0,end-start+1,4)
Mu[,1]=exp(model$summary.lincomb.derived$mean)
Mu[,2]=exp(model$summary.lincomb.derived$`0.025quant`)
Mu[,3]=exp(model$summary.lincomb.derived$`0.5quant`)
Mu[,4]=exp(model$summary.lincomb.derived$`0.975quant`)
#obtaining posterior R_t estimates
R[,1]=model$summary.fitted.values$mean
R[,2]=model$summary.fitted.values$`0.025quant`
R[,3]=model$summary.fitted.values$`0.5quant`
R[,4]=model$summary.fitted.values$`0.975quant`



#include imported case version
Imp_Infect=rep(0,end-start+1)
for(i in start:end)
{
  if(i>50){
    Imp_Infect[i-start+1]=sum(rev(w)*I_imp[(i-50):(i-1)])
  }else{
    Imp_Infect[i-start+1]=sum(rev(w[1:(i-1)])*I_imp[1:(i-1)])
  }
}
Imp_Infect=(Imp_Infect-mean(Imp_Infect))/sd(Imp_Infect)
data=data.frame(cbind(y,Infect,X_sd[start:end,],Imp_Infect))
data$idx = 1:nrow(data)
f1 <- as.formula(paste0("y ~ 1 + f(idx, model='iid') + ", paste0(colnames(X_sd), collapse = " + "),'+ Imp_Infect'))
ll=list()
for(i in 1:ncol(X_sd))
{
  ll[[colnames(X_sd)[i]]]=X_sd[start:end,i]
}
ll[['Imp_Infect']]=Imp_Infect
ll[["(Intercept)"]]=rep(1,end-start+1)
lc=inla.make.lincombs(ll)
model=inla(f1,data=data,family = "poisson",control.family=list(link='log'), E=Infect, lincomb =lc,control.compute=list(dic=TRUE, cpo=TRUE, return.marginals=TRUE, return.marginals.predictor=TRUE, config = TRUE),control.predictor = list(compute=TRUE))



#plot cases
library(grid)
library(scales)
plot_cases=function(region){
  start1=322
  xtk=cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31));xlb=c('J','F','M','A','M','J','J','A','S','O','N','D')
  out=which(xtk>=(end-start1+1))[-1]
  if (length(out)>0) {xtk=xtk[-out]; xlb=xlb[-(out-1)]}
  out=which.max(xtk>=(start-start1+1))
  if (out>2) {xtk=xtk[(out-1):length(xtk)]; xlb=xlb[(out-1):length(xlb)]}
  xtk[1]=start-start1
  end_date1=max(xtk)+start1-1
  pushViewport(plotViewport(c(4,4,1,1),xscale=c(xtk[1]+start1-1,end_date1),yscale=c(0,ym1+ym2))) #build another viewport in the first part
  grid.rect()
  grid.xaxis(at=xtk+start1-1,label=FALSE)
  grid.text(xlb,x=unit(0.5*(xtk[-1]+xtk[-length(xtk)])+start1-1,'native'),y=unit(-1,'lines'))
  grid.text('Month',y=unit(-2.5,'lines'))
  #grid.yaxis(at=ytk,gp=gpar(col='indianred', fontsize=7), main=FALSE)
  #grid.yaxis(at=ytk,label=ytk/r, gp=gpar(col='darkslateblue', fontsize=7))
  grid.yaxis(at=unit(ytk,'npc'),ylabel)
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  grid.lines(c(0,1),c(0.25,0.25))
  for(t in start:end)
  {
    grid.rect(x=t,y=unit(0.25,'npc'), width=1, height=unit(I[t]/ym1*0.75,'npc'),default.units = 'native',just=c('right','bottom'), gp=gpar(fill=scales::alpha('indianred',0.6),col=scales::alpha('indianred',0.1)))
    grid.rect(x=t,y=0, width=1, height=unit(I_imp[t]/ym2*0.25,"npc"),default.units = 'native',just=c('right','bottom'), gp=gpar(fill=scales::alpha('darkslateblue',0.6),col=scales::alpha('darkslateblue',0.1)))
  }
  if(region=='SG'){
    grid.text('Imported',x=unit(-3.5,'lines'),y=0.125,rot=90)
    grid.text('Autochthonous cases',x=unit(-3.5,'lines'),y=0.625,rot=90)
  }
  if(region=='NZ'){
    grid.text('Imported',x=unit(-3,'lines'),y=0.125,rot=90)
    grid.text('Autochthonous cases',x=unit(-3,'lines'),y=0.625,rot=90)
  }
  popViewport()
}

#SG until Sept 30, 2021
ym1=2100;ym2=60
ytk=c(seq(0,ym2-1,30)/ym2*0.25,0.25+seq(0,ym1,500)/ym1*0.75)
ylabel=c(seq(0,ym2-1,30),seq(0,ym1,500))
#SG
ym1=4800;ym2=200
ytk=c(seq(0,ym2-1,100)/ym2*0.25,0.25+seq(0,ym1,1000)/ym1*0.75)
ylabel=c(seq(0,ym2-1,100),seq(0,ym1,1000))
#NZ
ym1=240;ym2=30
ytk=c(seq(0,ym2-1,15)/ym2*0.25,0.25+seq(0,ym1,50)/ym1*0.75)
ylabel=c(seq(0,ym2-1,15),seq(0,ym1,50))
#NZ from August 19,2021
ym1=240;ym2=20
ytk=c(seq(0,ym2-1,10)/ym2*0.25,0.25+seq(0,ym1,50)/ym1*0.75)
ylabel=c(seq(0,ym2-1,10),seq(0,ym1,50))

png('cases.png',height=8,width=8,units='cm', res=300, pointsize=10) 
plot_cases(region)
dev.off()

